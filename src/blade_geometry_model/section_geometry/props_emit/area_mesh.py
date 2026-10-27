"""Thin adapter: section SDF components → shared 2-D quadtree area mesh."""
from __future__ import annotations
from typing import Any
import numpy as np
from blade_geometry_model.area_2d import AreaMesh2D, QuadtreeSpec, build_area_meshes_2d
from blade_geometry_model.area_2d.quadtree import (
    QuadTreeMesh,
    SupportBox,
    build_shared_quadtree,
    leaf_bounds,
    summarize_quadtree,
)
from blade_geometry_model.area_2d.quadrature import (
    _attach_center_phi,
    section_properties_all_labels,
)
from blade_geometry_model.area_2d.types import PhiSampler
from blade_geometry_model.section_geometry.midlines.section_midline_host import (
    resolve_section_host,
)
from blade_geometry_model.section_geometry.sdf.ir import expr_bbox
from blade_geometry_model.section_geometry.sdf.profiling import get_stats
from .section_evaluator import SectionFieldEvaluator

def _label_key(label: Any) -> str:
    val = getattr(label, "value", label)
    return str(val)

def _resolve_section_host(section: Any) -> Any:
    return resolve_section_host(section)

def _section_evaluator(section: Any, twist_rad: float = 0.0) -> SectionFieldEvaluator:
    del twist_rad
    return SectionFieldEvaluator.from_section(section)

def _section_label_keys(section: Any) -> list[str]:
    host = _resolve_section_host(section)
    raw_labels = list(getattr(host, "labels", ()) or ())
    if not raw_labels:
        raw_labels = list(host)
    return [_label_key(lbl) for lbl in raw_labels]

def section_phi_samplers(section: Any, twist_rad: float = 0.0) -> dict[str, PhiSampler]:
    """Build φ samplers for each labelled SDF component on ``section``."""

    evaluator = _section_evaluator(section, twist_rad)
    return evaluator.make_samplers(_section_label_keys(section))

def _clamp_support(box: SupportBox, margin_f: float, spec: QuadtreeSpec) -> SupportBox:
    x0, x1, y0, y1 = box
    return (
        max(spec.x_min, x0 - margin_f),
        min(spec.x_max, x1 + margin_f),
        max(spec.y_min, y0 - margin_f),
        min(spec.y_max, y1 + margin_f),
    )

def _component_supports(
    samplers: dict[str, PhiSampler],
    spec: QuadtreeSpec,
    *,
    exprs: dict[str, Any] | None = None,
    n: int = 64,
    margin: float | None = None,
) -> dict[str, SupportBox]:
    """Conservative support bbox per component.
    Prefers analytic IR-derived bounds (``expr_bbox``) — zero SDF evaluations —
    and falls back to a coarse uniform grid sample only for labels whose
    expression is unbounded (``HalfPlane``) or opaque (``CallableField``), or
    when no expression is available.
    """

    margin_f = float(margin if margin is not None else 3.0 * spec.min_cell_size)
    domain: SupportBox = (spec.x_min, spec.x_max, spec.y_min, spec.y_max)
    supports: dict[str, SupportBox] = {}
    fallback_labels: list[str] = []
    for lbl in samplers:
        expr = exprs.get(lbl) if exprs else None
        box = expr_bbox(expr) if expr is not None else None
        if box is None:
            fallback_labels.append(lbl)
            continue
        supports[lbl] = _clamp_support(box, margin_f, spec)
    if fallback_labels:
        xs = np.linspace(spec.x_min, spec.x_max, n)
        ys = np.linspace(spec.y_min, spec.y_max, n)
        xx, yy = np.meshgrid(xs, ys, indexing="xy")
        pts = np.column_stack([xx.ravel(), yy.ravel()])
        for lbl in fallback_labels:
            phi = np.asarray(samplers[lbl](pts), dtype=float).reshape(n, n)
            inside = phi < margin_f
            if not np.any(inside):
                supports[lbl] = domain
                continue
            ys_idx, xs_idx = np.where(inside)
            supports[lbl] = _clamp_support(
                (
                    float(xs[xs_idx].min()),
                    float(xs[xs_idx].max()),
                    float(ys[ys_idx].min()),
                    float(ys[ys_idx].max()),
                ),
                margin_f,
                spec,
            )
    return supports

def build_shared_section_quadtree(
    section: Any,
    spec: QuadtreeSpec,
    twist_rad: float = 0.0,
) -> QuadTreeMesh:
    """One quadtree refined to the union of all subcomponent boundaries."""

    evaluator = _section_evaluator(section, twist_rad)
    labels = _section_label_keys(section)
    samplers = evaluator.make_samplers(labels)
    supports = _component_supports(samplers, spec, exprs=evaluator.exprs_for(labels))
    multi_sampler = evaluator.make_multi_sampler(labels)
    mesh = build_shared_quadtree(
        spec, samplers, supports=supports, multi_sampler=multi_sampler
    )
    stats = get_stats()
    if stats is not None:
        stats.quadtree_leaf_count = int(leaf_bounds(mesh).shape[0])
    return mesh

def section_properties_via_quadtree(
    section: Any,
    spec: QuadtreeSpec,
    twist_rad: float = 0.0,
    labels: list[str] | None = None,
) -> tuple[dict[str, dict[str, float]], QuadTreeMesh, dict[str, object]]:
    """Compute component section properties on a shared adaptive quadtree.
    Parameters
    ----------
    labels :
        Optional subset of section labels. When provided, only those labels are
        refined / integrated (used by the hybrid residual-φ path).
    """

    evaluator = _section_evaluator(section, twist_rad)
    labels = list(labels) if labels is not None else _section_label_keys(section)
    if not labels:
        return {}, None, {"leaf_count": 0, "min_cell_width": None}  # type: ignore[return-value]
    samplers = evaluator.make_samplers(labels)
    supports = _component_supports(samplers, spec, exprs=evaluator.exprs_for(labels))
    multi_sampler = evaluator.make_multi_sampler(labels)
    mesh = build_shared_quadtree(
        spec, samplers, supports=supports, multi_sampler=multi_sampler
    )
    stats = get_stats()
    if stats is not None:
        stats.quadtree_leaf_count = int(leaf_bounds(mesh).shape[0])
    phi_by_label = getattr(mesh, "phi_by_label", None)
    if phi_by_label is None:
        raise RuntimeError("Shared quadtree missing phi_by_label after refinement.")
    phi_by_label = _attach_center_phi(
        mesh,
        phi_by_label,
        samplers,
        supports=supports,
        multi_sampler=multi_sampler,
    )
    props = section_properties_all_labels(mesh, phi_by_label, labels, samplers=samplers)
    qt_stats = summarize_quadtree(mesh)
    return props, mesh, qt_stats

def build_section_area_meshes(
    section: Any,
    spec: QuadtreeSpec,
    twist_rad: float = 0.0,
) -> AreaMesh2D:
    """Per-component quadtree meshes (legacy helper; prefer shared mesh for properties)."""

    return build_area_meshes_2d(section_phi_samplers(section, twist_rad), spec)
