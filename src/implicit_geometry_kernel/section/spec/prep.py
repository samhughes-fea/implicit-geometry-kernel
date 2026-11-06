"""Prep helpers and sized-section assembly for the geometry worker."""
from __future__ import annotations

from collections.abc import Collection, Mapping
from typing import Any, Sequence

from implicit_geometry_kernel.section.spec.layout_resolution import (
    _as_section_layout,
    _layout_web_fracs,
)
from implicit_geometry_kernel.section.spec.results import (
    ALL_OUTPUT_GATES,
    LayoutInput,
    normalize_output_gates,
)
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
from implicit_geometry_kernel.section.midlines.host import (
    SectionMidlineHost,
)
from implicit_geometry_kernel.sdf.shapes.airfoil import (
    AirfoilSDF,
)


def _compat_kwargs_from_zone_gauges(zone_gauges_m: Mapping[str, float]) -> dict[str, Any]:
    """One-way projection: group zone-id thicknesses into geometry wall gauges.

    Geometry stays gauge-agnostic and does not import ``section_gauge``.
    Callers that already hold a realization should pass ``zone_gauges_m``
    rather than collapsing to component maxima first.
    """
    items = {str(k): float(v) for k, v in dict(zone_gauges_m).items()}
    skins = [v for k, v in items.items() if k.startswith("skin:")]
    caps = [v for k, v in items.items() if k.startswith("cap:") or k.startswith("spar:")]
    webs = sorted(
        ((k, v) for k, v in items.items() if k.startswith("web:")),
        key=lambda kv: kv[0],
    )
    ss = float(items.get("skin:suction", max(skins) if skins else 0.0))
    ps = float(items.get("skin:pressure", ss))
    web_vals = [float(v) for _, v in webs]
    skin_kw: float | tuple[float, float] = ss if ss == ps else (ss, ps)
    if not web_vals:
        web_kw: float | list[float] = 0.0
    elif all(v == web_vals[0] for v in web_vals):
        web_kw = web_vals[0]
    else:
        web_kw = web_vals
    return {
        "skin_thickness": skin_kw,
        "spar_cap_thickness": float(max(caps) if caps else 0.0),
        "shear_web_thickness": web_kw,
    }


def _thickness_cache_token(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return tuple(round(float(x), 12) for x in value)
    return round(float(value), 12)


def section_geometry_prep_key(
    layout: LayoutInput,
    airfoil: AirfoilSDF,
    *,
    skin_thickness: float | None = None,
    spar_cap_thickness: float | None = None,
    shear_web_thickness: float | None = None,
    twist_rad: float = 0.0,
    props_backend: str | None = None,
    outputs: Collection[str] = ALL_OUTPUT_GATES,
    web_x_positions_m: Sequence[float] | None = None,
    zone_gauges_m: Mapping[str, float] | None = None,
    placement_digest_by_kind: Mapping[str, str] | None = None,
) -> tuple[Any, ...]:
    """Stable cache key for a geometry worker call (layout + airfoil + DVs + gates).

    ``web_x_positions_m`` participates: anchored positions change the built
    section without changing the layout key or any thickness, so omitting it
    would let two different sections collide on one cache entry.
    ``placement_digest_by_kind`` participates so a surface-aware policy change
    invalidates geometry caches (WP3 D3.3 / WP5 D5.3).
    """
    from implicit_geometry_kernel.section.spec.kwarg_aliases import (
        resolve_shear_web_thickness,
        resolve_skin_gauges,
        resolve_spar_cap_thickness,
    )

    if zone_gauges_m:
        derived = _compat_kwargs_from_zone_gauges(zone_gauges_m)
        skin_thickness = derived["skin_thickness"]
        spar_cap_thickness = derived["spar_cap_thickness"]
        shear_web_thickness = derived["shear_web_thickness"]

    t_skin_gauges = resolve_skin_gauges(
        skin_thickness=skin_thickness, required=True
    )
    t_spar_v = resolve_spar_cap_thickness(
        spar_cap_thickness=spar_cap_thickness, required=True
    )
    t_web_v = resolve_shear_web_thickness(
        shear_web_thickness=shear_web_thickness, required=True
    )
    chord = float(airfoil.chord or 0.0)
    naca = getattr(airfoil, "naca_code", None) or getattr(airfoil, "code", None)
    airfoil_id: Any
    if naca is not None:
        airfoil_id = ("naca", str(naca), chord)
    else:
        verts = airfoil.vertices
        if verts is not None:
            import numpy as np

            arr = np.asarray(verts, dtype=float)
            airfoil_id = ("verts", chord, int(arr.shape[0]), float(arr[0, 0]), float(arr[0, 1]))
        else:
            airfoil_id = ("airfoil", chord, id(airfoil))

    # Always key on parametric SectionLayout (catalogue specs are presets).
    layout_sot = layout if isinstance(layout, SectionLayout) else _as_section_layout(layout)
    cap = layout_sot.cap_placement
    layout_key = (
        "section_layout",
        layout_sot.catalogue_key,
        _layout_web_fracs(layout_sot),
        None
        if layout_sot.web_stations_m is None
        else tuple(float(x) for x in layout_sot.web_stations_m),
        str(layout_sot.orientation),
        str(cap.kind),
        str(cap.anchor),
        None if cap.half_width is None else float(cap.half_width),
        str(layout_sot.structural_family),
    )
    out_key = tuple(sorted(normalize_output_gates(outputs)))
    web_key = (
        None if web_x_positions_m is None
        else tuple(round(float(x), 12) for x in web_x_positions_m)
    )
    zone_key = (
        None
        if not zone_gauges_m
        else tuple(sorted((str(k), round(float(v), 12)) for k, v in zone_gauges_m.items()))
    )
    placement_key = (
        None
        if not placement_digest_by_kind
        else tuple(sorted((str(k), str(v)) for k, v in placement_digest_by_kind.items()))
    )
    return (
        layout_key,
        airfoil_id,
        _thickness_cache_token(t_skin_gauges),
        _thickness_cache_token(t_spar_v),
        _thickness_cache_token(t_web_v),
        float(twist_rad),
        props_backend,
        out_key,
        web_key,
        zone_key,
        placement_key,
    )


def build_section_with_thicknesses(
    airfoil: AirfoilSDF,
    layout: LayoutInput,
    *,
    twist_rad: float | None = None,
    twist_angle_rad: float | None = None,
    twist_angle: float | None = None,
    skin_thickness: float | None = None,
    spar_cap_thickness: float | None = None,
    shear_web_thickness: float | None = None,
    web_thickness: float | None = None,
    cap_height: float | None = None,
    web_x_positions_m: Sequence[float] | None = None,
    web_y_coords: Sequence[tuple[float, float]] | None = None,
    stage_sink: dict[str, float] | None = None,
    spar_cap_ply_pitch_m: float | None = None,
) -> SectionMidlineHost:
    """Build sized implicit section geometry for a layout and free DVs.

    Accepts catalogue :class:`SystemLayoutSpec` or parametric
    :class:`SectionLayout`. Every catalogue key (including ``0A`` / ``0B``)
    builds a :class:`MultiCellSection`. Prefer long thickness / twist names;
    ``cap_height`` / ``web_thickness`` remain accepted aliases.
    """
    from implicit_geometry_kernel.section.spec.kwarg_aliases import (
        resolve_shear_web_thickness,
        resolve_skin_gauges,
        resolve_spar_cap_thickness,
        resolve_twist_rad,
    )
    from implicit_geometry_kernel.section.assembly.from_layout import (
        multicell_from_layout,
    )

    return multicell_from_layout(
        airfoil,
        layout,
        twist_rad=resolve_twist_rad(
            twist_rad=twist_rad,
            twist_angle_rad=twist_angle_rad,
            twist_angle=twist_angle,
        ),
        skin_thickness=resolve_skin_gauges(  # type: ignore
            skin_thickness=skin_thickness, required=True
        ),
        spar_cap_thickness=resolve_spar_cap_thickness(
            spar_cap_thickness=spar_cap_thickness,
            cap_height=cap_height,
            required=True,
        ),
        shear_web_thickness=resolve_shear_web_thickness(  # type: ignore
            shear_web_thickness=shear_web_thickness,
            web_thickness=web_thickness,
            required=True,
        ),
        web_x_positions_m=web_x_positions_m,
        web_y_coords=web_y_coords,
        stage_sink=stage_sink,
        spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
    )


