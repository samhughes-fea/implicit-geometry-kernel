"""Build spar-cap wall recipes and SDF projections from CapPlacementRule.

Assemble authors frozen :class:`CapWall` recipes (G-RECIPE SoT), then projects
to SDF IR. Web-subtract remains an SDF-only post-step and is not encoded on
the recipe parameters.
"""

from __future__ import annotations

from typing import Any


from implicit_geometry_kernel.contract.errors import (
    InfeasibleSectionGeometry,
)
from implicit_geometry_kernel.sdf.nodes import (
    union_expr,
)

from implicit_geometry_kernel.section.layout.cap_placement import (
    CapPlacementKind,
    CapPlacementRule,
)
from implicit_geometry_kernel.section.layout.section_axes import max_thickness_chord_x, pitch_axis_x_from_le
from implicit_geometry_kernel.section.assembly.cap_feasibility import ExtentAnchor, cap_feasibility_clamp
from implicit_geometry_kernel.section.assembly.build_spar_caps_pads import (
    _cap_wall,
    _label_str,
    _project_cap,
    _register_cap_side,
    _wrap_expr_cap,
    snap_cap_heights,
)
from implicit_geometry_kernel.section.assembly.clip_algebra import (
    _parallel_web_strip_clip,
    flapwise_strip_params,
    web_station_projection,
)
from implicit_geometry_kernel.section.assembly.labels import ComponentLabel
from implicit_geometry_kernel.section.assembly.web_layout import inner_y_at_x
from implicit_geometry_kernel.section.walls.walls import CapWall


def build_spar_caps(
    airfoil_sdf,
    skin_thickness,
    xs,
    anchors,
    cap_h_upper,
    cap_h_lower,
    web_components,
    *,
    cap_placement: CapPlacementRule,
    airfoil_expr_node,
    all_flapwise: bool,
    c_station,
    nx,
    ny,
    N: int,
    pitch_fr: float,
    fhw: float,
    dhw: float,
    y_search,
    phi_inner_ref,
    spar_cap_ply_pitch_m: float | None = None,
) -> tuple[Any | None, Any | None, dict[str, Any], dict[str, CapWall]]:
    """Return (cap_upper, cap_lower, unrotated_SparCap, CapWall recipes).

    CapWall recipes are authored first from :class:`CapPlacementRule`; SDF
    components are projections. Family letters are not accepted here.

    ``spar_cap_ply_pitch_m`` (optional): after cavity feasibility reduces a
    continuous height, floor onto this pitch so CapWall / midline thickness
    stays on the ply pitch grid. Geometry remains material-agnostic —
    callers pass the catalogue spar-cap ``PlyGrid.pitch`` when known.
    """
    rule = cap_placement
    build = _CAP_BUILDERS.get(rule.kind)
    if build is None:
        raise ValueError(f"Unknown CapPlacementRule.kind {rule.kind!r}")
    if rule.kind == "none":
        return None, None, {}, {}

    web_exprs = [w.expr for w in web_components]
    web_union = union_expr(web_exprs) if web_exprs else None
    unrotated: dict[str, Any] = {}
    recipes: dict[str, CapWall] = {}
    cap_upper, cap_lower = build(
        rule=rule,
        airfoil_sdf=airfoil_sdf,
        skin_thickness=skin_thickness,
        xs=xs,
        anchors=anchors,
        cap_h_upper=cap_h_upper,
        cap_h_lower=cap_h_lower,
        N=N,
        pitch_fr=pitch_fr,
        fhw=fhw,
        dhw=dhw,
        y_search=y_search,
        phi_inner_ref=phi_inner_ref,
        spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
        all_flapwise=all_flapwise,
        c_station=c_station,
        nx=nx,
        ny=ny,
        recipes=recipes,
        unrotated=unrotated,
        airfoil_expr_node=airfoil_expr_node,
        web_union=web_union,
    )
    return cap_upper, cap_lower, unrotated, recipes


def _caps_none(**_: Any) -> tuple[Any, Any]:
    return None, None


def _caps_continuous_box(
    *,
    xs,
    N: int,
    airfoil_sdf,
    skin_thickness,
    cap_h_upper,
    cap_h_lower,
    y_search,
    phi_inner_ref,
    spar_cap_ply_pitch_m,
    all_flapwise,
    c_station,
    nx,
    ny,
    recipes,
    unrotated,
    airfoil_expr_node,
    web_union,
    **_: Any,
) -> tuple[Any, Any]:
    if N < 2:
        raise InfeasibleSectionGeometry(
            f"continuous_box cap placement requires at least 2 webs (got N={N}). There is no 1D catalogue layout."
        )
    x_cap_start = xs[0]
    x_cap_end = xs[-1]
    x_cap_start, x_cap_end, cap_h_upper, cap_h_lower, _clamp_info = cap_feasibility_clamp(
        airfoil_sdf,
        skin_thickness,
        x_cap_start,
        x_cap_end,
        cap_h_upper,
        cap_h_lower,
        y_search=y_search,
        phi_inner_ref=phi_inner_ref,
        aft_half_depth_rule=True,
    )
    cap_h_upper, cap_h_lower = snap_cap_heights(cap_h_upper, cap_h_lower, spar_cap_ply_pitch_m)
    cap_strip = None
    cap_fw: dict[str, Any] = {}
    if all_flapwise and c_station is not None:
        c_lo_cap = float(c_station[0])
        c_hi_cap = float(c_station[-1])
        d_stat = abs(c_lo_cap - c_hi_cap)
        chord_scale = float(getattr(airfoil_sdf, "chord", 1.0))
        if d_stat > 1e-9 * max(chord_scale, 1e-9):
            cap_strip = _parallel_web_strip_clip(nx, ny, c_lo_cap, c_hi_cap)
            cap_fw = flapwise_strip_params(nx, ny, c_lo_cap, c_hi_cap)
    _camber = _clamp_info.get("camber_xy")
    wall_u = _cap_wall(
        skin_thickness=skin_thickness,
        x_start=x_cap_start,
        x_end=x_cap_end,
        cap_height=cap_h_upper,
        surface="suction_surface",
        fw=cap_fw,
        camber_xy=_camber,
    )
    wall_l = _cap_wall(
        skin_thickness=skin_thickness,
        x_start=x_cap_start,
        x_end=x_cap_end,
        cap_height=cap_h_lower,
        surface="pressure_surface",
        fw=cap_fw,
        camber_xy=_camber,
    )
    cap_upper = _register_cap_side(
        wall=wall_u,
        label=ComponentLabel.SPAR_CAP_SUCTION_SURFACE,
        recipes=recipes,
        unrotated=unrotated,
        airfoil_sdf=airfoil_sdf,
        airfoil_expr_node=airfoil_expr_node,
        camber_xy=_camber,
        strip_clip=cap_strip,
        web_union=web_union,
    )
    cap_lower = _register_cap_side(
        wall=wall_l,
        label=ComponentLabel.SPAR_CAP_PRESSURE_SURFACE,
        recipes=recipes,
        unrotated=unrotated,
        airfoil_sdf=airfoil_sdf,
        airfoil_expr_node=airfoil_expr_node,
        camber_xy=_camber,
        strip_clip=cap_strip,
        web_union=web_union,
    )
    return cap_upper, cap_lower


def _caps_fixed_band(
    *,
    rule,
    airfoil_sdf,
    skin_thickness,
    cap_h_upper,
    cap_h_lower,
    y_search,
    phi_inner_ref,
    spar_cap_ply_pitch_m,
    all_flapwise,
    c_station,
    nx,
    ny,
    pitch_fr,
    fhw,
    recipes,
    unrotated,
    airfoil_expr_node,
    web_union,
    **_: Any,
) -> tuple[Any, Any]:
    window = rule.chord_window
    extent_anchor: ExtentAnchor
    if window is not None:
        chord = float(getattr(airfoil_sdf, "chord", 1.0))
        lo_f, hi_f = float(window[0]), float(window[1])
        x_lo = min(lo_f, hi_f) * chord
        x_hi = max(lo_f, hi_f) * chord
        x_c = 0.5 * (x_lo + x_hi)
        extent_anchor = "deepest"
    else:
        anchor = rule.anchor or "max_thickness"
        half = float(rule.half_width) if rule.half_width is not None else float(fhw)
        if anchor == "pitching":
            x_c = pitch_axis_x_from_le(airfoil_sdf, pitch_fr)
        else:
            x_c = max_thickness_chord_x(airfoil_sdf)
        x_lo = x_c - half
        x_hi = x_c + half
        extent_anchor = "band_center"
    chord_m = float(getattr(airfoil_sdf, "chord", 1.0))
    if rule.side_windows is not None:
        # each side's cap in its own window (a published layup whose suction and
        # pressure caps differ); each band is clamped with both cap heights
        bands = [(min(w) * chord_m, max(w) * chord_m) for w in rule.side_windows]
        extent_anchor = "deepest"
    else:
        bands = [(x_lo, x_hi), (x_lo, x_hi)]
    sides: list[tuple] = []
    for b_lo, b_hi in bands:
        if sides and rule.side_windows is None:
            sides.append(sides[0])  # one band for both sides: clamped once
            break
        b_c = float(x_c) if rule.side_windows is None else 0.5 * (b_lo + b_hi)
        sides.append((b_c, *cap_feasibility_clamp(
            airfoil_sdf,
            skin_thickness,
            b_lo,
            b_hi,
            cap_h_upper,
            cap_h_lower,
            y_search=y_search,
            phi_inner_ref=phi_inner_ref,
            aft_half_depth_rule=False,
            extent_anchor=extent_anchor,
            center_x=b_c,
        )))
    (c_u, x_lo_u, x_hi_u, h_u, _, info_u), (c_l, x_lo_l, x_hi_l, _, h_l, info_l) = sides
    cap_h_upper, cap_h_lower = snap_cap_heights(h_u, h_l, spar_cap_ply_pitch_m)
    out = []
    for surface, label, centre, lo, hi, height, info in (
            ("suction_surface", ComponentLabel.SPAR_CAP_SUCTION_SURFACE, c_u, x_lo_u, x_hi_u,
             cap_h_upper, info_u),
            ("pressure_surface", ComponentLabel.SPAR_CAP_PRESSURE_SURFACE, c_l, x_lo_l, x_hi_l,
             cap_h_lower, info_l)):
        camber = info.get("camber_xy")
        strip = None
        strip_fw: dict[str, Any] = {}
        if all_flapwise and c_station is not None:
            y_tc, y_bc = inner_y_at_x(
                airfoil_sdf,
                skin_thickness,
                centre,
                y_search=y_search,
                phi_inner=phi_inner_ref,
            )
            cy_c = 0.5 * (y_tc + y_bc)
            c_lo = web_station_projection(nx, ny, lo, cy_c)
            c_hi = web_station_projection(nx, ny, hi, cy_c)
            strip = _parallel_web_strip_clip(nx, ny, c_lo, c_hi)
            strip_fw = flapwise_strip_params(nx, ny, c_lo, c_hi)
        wall = _cap_wall(
            skin_thickness=skin_thickness,
            x_start=lo,
            x_end=hi,
            cap_height=height,
            surface=surface,
            fw=strip_fw,
            camber_xy=camber,
        )
        out.append(_register_cap_side(
            wall=wall,
            label=label,
            recipes=recipes,
            unrotated=unrotated,
            airfoil_sdf=airfoil_sdf,
            airfoil_expr_node=airfoil_expr_node,
            camber_xy=camber,
            strip_clip=strip,
            web_union=web_union,
        ))
    return out[0], out[1]


def _caps_per_web_pads(
    *,
    airfoil_sdf,
    skin_thickness,
    xs,
    anchors,
    cap_h_upper,
    cap_h_lower,
    N: int,
    dhw,
    y_search,
    phi_inner_ref,
    spar_cap_ply_pitch_m,
    all_flapwise,
    c_station,
    nx,
    ny,
    recipes,
    unrotated,
    airfoil_expr_node,
    web_union,
    **_: Any,
) -> tuple[Any, Any]:
    from implicit_geometry_kernel.section.assembly.build_spar_caps_pads import fill_per_web_pads

    return fill_per_web_pads(
        airfoil_sdf=airfoil_sdf,
        skin_thickness=skin_thickness,
        xs=xs,
        anchors=anchors,
        cap_h_upper=cap_h_upper,
        cap_h_lower=cap_h_lower,
        N=N,
        dhw=dhw,
        y_search=y_search,
        phi_inner_ref=phi_inner_ref,
        spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
        all_flapwise=all_flapwise,
        c_station=c_station,
        nx=nx,
        ny=ny,
        recipes=recipes,
        unrotated=unrotated,
        airfoil_expr_node=airfoil_expr_node,
        web_union=web_union,
        _cap_wall=_cap_wall,
        _project_cap=_project_cap,
        _wrap_expr_cap=_wrap_expr_cap,
        _label_str=_label_str,
    )


_CAP_BUILDERS: dict[CapPlacementKind, Any] = {
    "none": _caps_none,
    "continuous_box": _caps_continuous_box,
    "fixed_band": _caps_fixed_band,
    "per_web_pads": _caps_per_web_pads,
}


__all__ = ["build_spar_caps"]
