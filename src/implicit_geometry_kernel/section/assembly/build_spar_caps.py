"""Build spar-cap wall recipes and SDF projections from CapPlacementRule.

Assemble authors frozen :class:`CapWall` recipes (G-RECIPE SoT), then projects
to SDF IR. Web-subtract remains an SDF-only post-step and is not encoded on
the recipe parameters.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from implicit_geometry_kernel.contract.errors import (
    InfeasibleSectionGeometry,
)
from implicit_geometry_kernel.section.walls.walls import (
    CapWall,
)
from implicit_geometry_kernel.sdf.nodes import (
    Expr,
    Subtract,
    union_expr,
)

from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
from implicit_geometry_kernel.section.layout.section_axes import max_thickness_chord_x, pitch_axis_x_from_le
from implicit_geometry_kernel.section.assembly.cap_feasibility import ExtentAnchor, cap_feasibility_clamp
from implicit_geometry_kernel.section.assembly.cap_ply_snap import snap_cap_heights
from implicit_geometry_kernel.section.assembly.clip_algebra import (
    _parallel_web_strip_clip,
    flapwise_strip_params,
    web_station_projection,
)
from implicit_geometry_kernel.section.walls.spar_cap import (
    SparCap,
)
from implicit_geometry_kernel.section.assembly.labels import ComponentLabel
from implicit_geometry_kernel.section.assembly.web_layout import inner_y_at_x


def _label_str(label: Any) -> str:
    return str(getattr(label, "value", label))


def _wrap_expr_cap(base_cap: SparCap, expr: Expr) -> SparCap:
    return base_cap.with_expr(expr)


def _camber_xy_tuple(camber_xy) -> tuple[tuple[float, float], ...] | None:
    """Freeze feasibility camber samples for CapWall (hashable recipe field)."""
    if camber_xy is None:
        return None
    arr = np.asarray(camber_xy, dtype=float)
    if arr.ndim != 2 or arr.shape[0] < 2 or arr.shape[1] < 2:
        return None
    return tuple((float(x), float(y)) for x, y in arr[:, :2])


def _cap_wall(
    *,
    skin_thickness: float,
    x_start: float,
    x_end: float,
    cap_height: float,
    surface: str,
    twist_angle: float = 0.0,
    fw: dict[str, Any] | None = None,
    camber_xy=None,
) -> CapWall | None:
    """Build a CapWall recipe, or ``None`` when snapped height is geometric run-out."""
    if float(cap_height) <= 0.0:
        return None
    fw = fw or {}
    return CapWall(
        skin_thickness=float(skin_thickness),
        x_start=float(x_start),
        x_end=float(x_end),
        cap_height=float(cap_height),
        surface=str(surface),
        twist_angle=float(twist_angle),
        flapwise_nx=fw.get("flapwise_nx"),
        flapwise_ny=fw.get("flapwise_ny"),
        flapwise_c_lo=fw.get("flapwise_c_lo"),
        flapwise_c_hi=fw.get("flapwise_c_hi"),
        flapwise_clip_mode=fw.get("flapwise_clip_mode"),
        camber_xy=_camber_xy_tuple(camber_xy),
    )


def _project_cap(
    wall: CapWall,
    airfoil_sdf,
    *,
    airfoil_expr_node,
    camber_xy=None,
    strip_clip=None,
) -> SparCap:
    cap = wall.to_sdf_component(
        airfoil_sdf, camber_xy=camber_xy, airfoil_expr_node=airfoil_expr_node
    )
    if strip_clip is not None:
        cap._strip_clip = strip_clip
    return cap


def _register_cap_side(
    *,
    wall: CapWall | None,
    label: ComponentLabel,
    recipes: dict[str, Any],
    unrotated: dict[Any, Any],
    airfoil_sdf,
    airfoil_expr_node,
    camber_xy,
    strip_clip,
    web_union,
) -> SparCap | None:
    """Project and register one cap side; skip geometric run-out (``wall is None``)."""
    if wall is None:
        return None
    recipes[_label_str(label)] = wall
    cap = _project_cap(
        wall,
        airfoil_sdf,
        airfoil_expr_node=airfoil_expr_node,
        camber_xy=camber_xy,
        strip_clip=strip_clip,
    )
    if web_union is not None:
        cap = _wrap_expr_cap(cap, Subtract(cap.expr, web_union))
        cap.footprint_expr = Subtract(cap.footprint_expr, web_union)
    unrotated[label] = cap
    return cap


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

    if rule.kind == "none":
        return None, None, {}, {}

    web_exprs = [w.expr for w in web_components]
    web_union = union_expr(web_exprs) if web_exprs else None
    unrotated: dict[str, Any] = {}
    recipes: dict[str, CapWall] = {}
    cap_upper = cap_lower = None

    if rule.kind == "continuous_box":
        # No 1D catalogue layout: a continuous box needs fore and aft webs.
        # Resolve must step D -> 1C (per_web_pads) before assemble.
        if N < 2:
            raise InfeasibleSectionGeometry(
                "continuous_box cap placement requires at least 2 webs "
                f"(got N={N}). There is no 1D catalogue layout."
            )
        x_cap_start = xs[0]
        x_cap_end = xs[-1]
        x_cap_start, x_cap_end, cap_h_upper, cap_h_lower, _clamp_info = (
            cap_feasibility_clamp(
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
        )
        cap_h_upper, cap_h_lower = snap_cap_heights(
            cap_h_upper, cap_h_lower, spar_cap_ply_pitch_m
        )
        cap_strip = None
        cap_fw: dict[str, Any] = {}
        if all_flapwise and c_station is not None:
            # Web midplane stations are the strip SoT (same as sandwich cores).
            # Do not recompute from camber at x_cap — that drifts under twist.
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
    elif rule.kind == "fixed_band":
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
        x_lo, x_hi, cap_h_upper, cap_h_lower, _clamp_b = cap_feasibility_clamp(
            airfoil_sdf,
            skin_thickness,
            x_lo,
            x_hi,
            cap_h_upper,
            cap_h_lower,
            y_search=y_search,
            phi_inner_ref=phi_inner_ref,
            aft_half_depth_rule=False,
            extent_anchor=extent_anchor,
            center_x=float(x_c),
        )
        cap_h_upper, cap_h_lower = snap_cap_heights(
            cap_h_upper, cap_h_lower, spar_cap_ply_pitch_m
        )
        _camber_b = _clamp_b.get("camber_xy")
        strip_b = None
        strip_b_fw: dict[str, Any] = {}
        if all_flapwise and c_station is not None:
            y_tc, y_bc = inner_y_at_x(
                airfoil_sdf,
                skin_thickness,
                x_c,
                y_search=y_search,
                phi_inner=phi_inner_ref,
            )
            cy_c = 0.5 * (y_tc + y_bc)
            c_lo = web_station_projection(nx, ny, x_lo, cy_c)
            c_hi = web_station_projection(nx, ny, x_hi, cy_c)
            strip_b = _parallel_web_strip_clip(nx, ny, c_lo, c_hi)
            strip_b_fw = flapwise_strip_params(nx, ny, c_lo, c_hi)
        wall_u = _cap_wall(
            skin_thickness=skin_thickness,
            x_start=x_lo,
            x_end=x_hi,
            cap_height=cap_h_upper,
            surface="suction_surface",
            fw=strip_b_fw,
            camber_xy=_camber_b,
        )
        wall_l = _cap_wall(
            skin_thickness=skin_thickness,
            x_start=x_lo,
            x_end=x_hi,
            cap_height=cap_h_lower,
            surface="pressure_surface",
            fw=strip_b_fw,
            camber_xy=_camber_b,
        )
        cap_upper = _register_cap_side(
            wall=wall_u,
            label=ComponentLabel.SPAR_CAP_SUCTION_SURFACE,
            recipes=recipes,
            unrotated=unrotated,
            airfoil_sdf=airfoil_sdf,
            airfoil_expr_node=airfoil_expr_node,
            camber_xy=_camber_b,
            strip_clip=strip_b,
            web_union=web_union,
        )
        cap_lower = _register_cap_side(
            wall=wall_l,
            label=ComponentLabel.SPAR_CAP_PRESSURE_SURFACE,
            recipes=recipes,
            unrotated=unrotated,
            airfoil_sdf=airfoil_sdf,
            airfoil_expr_node=airfoil_expr_node,
            camber_xy=_camber_b,
            strip_clip=strip_b,
            web_union=web_union,
        )
    elif rule.kind == "per_web_pads":
        from implicit_geometry_kernel.section.assembly.build_spar_caps_pads import fill_per_web_pads

        cap_upper, cap_lower = fill_per_web_pads(
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

    return cap_upper, cap_lower, unrotated, recipes


__all__ = ["build_spar_caps"]
