"""Section property backend wiring for SectionPropertiesReport."""

from __future__ import annotations

import time
from typing import Any, Literal

import numpy as np

from blade_geometry_model.section_geometry.props_emit.area_mesh import (
    section_properties_via_quadtree,
)
from blade_geometry_model.section_geometry.constructive import (
    classify_property_labels,
    emit_constructive_polygons,
    section_properties_from_polygons,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    PropsBackend,
    props_backend_from_env,
    quadtree_spec_for_airfoil,
    quadtree_spec_to_dict,
)
from blade_geometry_model.area_2d.spec import QuadtreeSpec
from blade_geometry_model.section_geometry.sdf.cache import SectionEvalCache
from blade_geometry_model.section_geometry.sdf.profiling import get_stats, timed_block


PropsBackendArg = Literal["grid", "quadtree", "hybrid"] | PropsBackend


def _label_key(label: Any) -> str:
    val = getattr(label, "value", label)
    return str(val)


def _section_labels(section_geometry: Any) -> list[str]:
    raw = list(getattr(section_geometry, "labels", ()) or ())
    if not raw:
        raw = list(section_geometry)
    return [_label_key(lbl) for lbl in raw]


def _airfoil_bbox_fallback(grid) -> Any:
    xmin, xmax = float(grid.X.min()), float(grid.X.max())
    ymin, ymax = float(grid.Y.min()), float(grid.Y.max())
    return type(
        "_BBox",
        (),
        {
            "vertices": np.array(
                [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]],
                dtype=float,
            )
        },
    )()


def _compute_hybrid(
    section_geometry: Any,
    *,
    spec: QuadtreeSpec,
    twist_rad: float,
    labels: list[str],
) -> tuple[dict[str, dict[str, float]], Any, dict[str, Any]]:
    constructive_labels, residual_labels = classify_property_labels(labels)
    meta_extra: dict[str, Any] = {
        "constructive_labels": list(constructive_labels),
        "wall_solid_labels": [k for k in constructive_labels if not k.startswith("core")],
        "residual_labels": list(residual_labels),
        "residual_maturity": {
            "core": "constructive polygon path (cavity - wall solids + bay clip)",
            "leading_edge_insert": "constructive polygon path (default)",
            "trailing_edge_insert": "constructive polygon path (default)",
        },
    }

    emit_timings: dict[str, float] = {}
    t0 = time.perf_counter()
    polys = emit_constructive_polygons(section_geometry, stage_sink=emit_timings)
    t_emit = time.perf_counter()
    polys = {
        k: v
        for k, v in polys.items()
        if k in constructive_labels
        and v is not None
        and not getattr(v, "is_empty", True)
        and float(getattr(v, "area", 0.0) or 0.0) > 1e-16
    }
    t_filter = time.perf_counter()
    constructive_props = section_properties_from_polygons(polys)
    t_shoelace = time.perf_counter()
    constructive_s = t_shoelace - t0
    meta_extra["constructive_s"] = float(constructive_s)
    meta_extra["properties_emit"] = float(t_emit - t0)
    meta_extra["properties_filter"] = float(t_filter - t_emit)
    meta_extra["properties_shoelace"] = float(t_shoelace - t_filter)
    meta_extra["polygon_count"] = int(len(polys))
    for k, v in emit_timings.items():
        meta_extra[str(k)] = float(v)

    mesh = None
    qt_stats: dict[str, object] = {}
    residual_props: dict[str, dict[str, float]] = {}
    t1 = time.perf_counter()
    if residual_labels:
        residual_props, mesh, qt_stats = section_properties_via_quadtree(
            section_geometry,
            spec,
            twist_rad=float(twist_rad),
            labels=list(residual_labels),
        )
    residual_s = time.perf_counter() - t1
    meta_extra["residual_quadtree_s"] = float(residual_s)
    meta_extra["quadtree_stats"] = qt_stats

    props_map = dict(constructive_props)
    props_map.update(residual_props)
    # Ensure every section label has an entry.
    for label in labels:
        if label not in props_map:
            props_map[label] = {
                "area": 0.0,
                "cx": float("nan"),
                "cy": float("nan"),
                "Ixx": float("nan"),
                "Iyy": float("nan"),
                "Ixy": float("nan"),
                "r_gyr_x": float("nan"),
                "r_gyr_y": float("nan"),
                "label": label,
            }
        else:
            props_map[label]["label"] = label
    meta_extra["constructive_polygons"] = {
        k: "emitted" for k in polys
    }
    return props_map, mesh, meta_extra


def compute_section_properties(
    section_geometry: Any,
    *,
    backend: PropsBackendArg | None = None,
    grid=None,
    eval_cache: SectionEvalCache | None = None,
    quadtree_spec: QuadtreeSpec | None = None,
    chord: float = 1.0,
    twist_rad: float = 0.0,
    airfoil_for_spec: Any | None = None,
    with_mesh: bool = False,
) -> tuple[dict[str, dict[str, float]], dict[str, Any]]:
    """Compute section properties for all labels using grid, quadtree, or hybrid.

    When ``with_mesh=True`` the return is a 3-tuple ``(data, meta, mesh)`` where
    ``mesh`` is the shared quadtree (full-label for ``quadtree``, residual-only
    for ``hybrid``) and ``None`` for the grid backend.
    """
    backend_use: PropsBackend = backend or props_backend_from_env()  # type: ignore[assignment]
    meta: dict[str, Any] = {"backend": str(backend_use)}
    labels = _section_labels(section_geometry)

    if backend_use == "grid":
        if grid is None:
            raise ValueError("grid is required when backend='grid'.")
        with timed_block("properties"):
            data: dict[str, dict[str, float]] = {}
            for label in labels:
                if eval_cache is None:
                    phi = grid.eval(section_geometry[label])
                else:
                    phi = eval_cache.get_or_eval_with_owner(
                        label, section_geometry[label], grid, owner=section_geometry
                    )
                props = grid.section_properties(phi)
                props["label"] = label
                data[label] = props
        meta["grid"] = {"nx": int(grid.nx), "ny": int(grid.ny)}
        return (data, meta, None) if with_mesh else (data, meta)

    af = airfoil_for_spec
    if af is None and grid is not None:
        af = _airfoil_bbox_fallback(grid)
    spec = quadtree_spec or quadtree_spec_for_airfoil(af, chord=float(chord))
    meta["quadtree_spec"] = quadtree_spec_to_dict(spec)

    if backend_use == "hybrid":
        with timed_block("properties"):
            props_map, mesh, hybrid_meta = _compute_hybrid(
                section_geometry,
                spec=spec,
                twist_rad=float(twist_rad),
                labels=labels,
            )
        meta.update(hybrid_meta)
        stats = get_stats()
        if stats is not None:
            meta["profiling"] = stats.to_dict()
        return (props_map, meta, mesh) if with_mesh else (props_map, meta)

    with timed_block("properties"):
        props_map, mesh, qt_stats = section_properties_via_quadtree(
            section_geometry, spec, twist_rad=float(twist_rad)
        )
    for label in labels:
        props_map[label]["label"] = label
    meta["quadtree_stats"] = qt_stats
    stats = get_stats()
    if stats is not None:
        meta["profiling"] = stats.to_dict()
    return (props_map, meta, mesh) if with_mesh else (props_map, meta)
