"""Build spar-cap wall recipes and SDF projections from CapPlacementRule.

Assemble authors frozen :class:`CapWall` recipes (G-RECIPE SoT), then projects
to SDF IR. Web-subtract remains an SDF-only post-step and is not encoded on
the recipe parameters.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from blade_geometry_model.section_contract.errors import (
    InfeasibleSectionGeometry,
)
from blade_geometry_model.section_geometry.recipes.walls import (
    CapWall,
)
from blade_geometry_model.section_geometry.sdf.ir import (
    Expr,
    Subtract,
    union_expr,
)

from ..layout.cap_placement import CapPlacementRule
from ..layout.section_axes import max_thickness_chord_x, pitch_axis_x_from_le
from .cap_feasibility import cap_feasibility_clamp
from .cap_ply_snap import snap_cap_heights
from .clip_algebra import _parallel_web_strip_clip, web_station_projection
from blade_geometry_model.section_geometry.recipes.sdf_ir.spar_cap import (
    SparCap,
)
from .labels import ComponentLabel
from .web_layout import inner_y_at_x


def _label_str(label: Any) -> str:
    return str(getattr(label, "value", label))


def _flapwise_strip_params(nx, ny, c_lo, c_hi) -> dict[str, Any]:
    """Kwargs mirroring a ``_parallel_web_strip_clip`` for the Shapely emitter."""
    return {
        "flapwise_nx": float(nx),
        "flapwise_ny": float(ny),
        "flapwise_c_lo": float(min(c_lo, c_hi)),
        "flapwise_c_hi": float(max(c_lo, c_hi)),
        "flapwise_clip_mode": "strip",
    }


def _wrap_expr_cap(base_cap: SparCap, expr: Expr) -> SparCap:
    out = SparCap.__new__(SparCap)
    out.surface = base_cap.surface
    out.x_start = base_cap.x_start
    out.x_end = base_cap.x_end
    out.cap_height = base_cap.cap_height
    out.skin_thickness = base_cap.skin_thickness
    out.twist_angle = base_cap.twist_angle
    out._af = base_cap._af
    out._strip_clip = base_cap._strip_clip
    out._camber_xy = getattr(base_cap, "_camber_xy", None)
    out._camber_surface = getattr(base_cap, "_camber_surface", None)
    for _fw in (
        "flapwise_nx",
        "flapwise_ny",
        "flapwise_c_lo",
        "flapwise_c_hi",
        "flapwise_clip_mode",
    ):
        setattr(out, _fw, getattr(base_cap, _fw, None))
    out.expr = expr
    out._expr = expr
    return out


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
    stays on the ply pitch grid (BL-019). Geometry remains material-agnostic —
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
            # Do not recompute from camber at x_cap — that drifts under twist (SG-011).
            c_lo_cap = float(c_station[0])
            c_hi_cap = float(c_station[-1])
            d_stat = abs(c_lo_cap - c_hi_cap)
            chord_scale = float(getattr(airfoil_sdf, "chord", 1.0))
            if d_stat > 1e-9 * max(chord_scale, 1e-9):
                cap_strip = _parallel_web_strip_clip(nx, ny, c_lo_cap, c_hi_cap)
                cap_fw = _flapwise_strip_params(nx, ny, c_lo_cap, c_hi_cap)
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
        anchor = rule.anchor or "max_thickness"
        half = float(rule.half_width) if rule.half_width is not None else float(fhw)
        if anchor == "pitching":
            x_c = pitch_axis_x_from_le(airfoil_sdf, pitch_fr)
        else:
            x_c = max_thickness_chord_x(airfoil_sdf)
        x_lo = x_c - half
        x_hi = x_c + half
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
            extent_anchor="band_center",
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
            strip_b_fw = _flapwise_strip_params(nx, ny, c_lo, c_hi)
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
        if N < 1:
            raise ValueError(
                "per_web_pads cap placement requires at least one web "
                f"(got N={N})."
            )
        # First pass: clamp each pad at the requested height to learn the
        # shared (min) feasible height across webs that can carry a centred pad.
        _clamped_c = [
            cap_feasibility_clamp(
                airfoil_sdf,
                skin_thickness,
                xs[i] - dhw,
                xs[i] + dhw,
                cap_h_upper,
                cap_h_lower,
                y_search=y_search,
                phi_inner_ref=phi_inner_ref,
                aft_half_depth_rule=False,
                extent_anchor="band_center",
                center_x=float(xs[i]),
            )
            for i in range(N)
        ]
        _live = [
            c
            for c in _clamped_c
            if float(c[2]) > 0.0
            and float(c[3]) > 0.0
            and (float(c[1]) - float(c[0])) > 1e-12
            and not bool(c[4].get("pad_run_out"))
        ]
        if _live:
            cap_h_upper = min(float(c[2]) for c in _live)
            cap_h_lower = min(float(c[3]) for c in _live)
        elif _clamped_c:
            # Every pad ran out — keep the clamp heights (all zero) so CapWall
            # recipes are omitted rather than inventing a gauge.
            cap_h_upper = min(float(c[2]) for c in _clamped_c)
            cap_h_lower = min(float(c[3]) for c in _clamped_c)
        cap_h_upper, cap_h_lower = snap_cap_heights(
            cap_h_upper, cap_h_lower, spar_cap_ply_pitch_m
        )
        # Second pass (SG-013): re-clamp extents at the shared height so pads
        # that were over-trimmed at the taller request can grow back
        # symmetrically within +/- dhw about each web.
        _clamped_c = [
            cap_feasibility_clamp(
                airfoil_sdf,
                skin_thickness,
                xs[i] - dhw,
                xs[i] + dhw,
                cap_h_upper,
                cap_h_lower,
                y_search=y_search,
                phi_inner_ref=phi_inner_ref,
                aft_half_depth_rule=False,
                extent_anchor="band_center",
                center_x=float(xs[i]),
            )
            for i in range(N)
        ]

        caps_u = []
        caps_l = []
        for i in range(N):
            x_lo, x_hi = float(_clamped_c[i][0]), float(_clamped_c[i][1])
            if (x_hi - x_lo) < 1e-12 or bool(_clamped_c[i][4].get("pad_run_out")):
                continue
            _camber_c = _clamped_c[i][4].get("camber_xy")
            strip_c = None
            strip_c_fw: dict[str, Any] = {}
            if all_flapwise and c_station is not None:
                cy_i = 0.5 * (anchors[i][0] + anchors[i][1])
                c_lo = web_station_projection(nx, ny, x_lo, cy_i)
                c_hi = web_station_projection(nx, ny, x_hi, cy_i)
                strip_c = _parallel_web_strip_clip(nx, ny, c_lo, c_hi)
                strip_c_fw = _flapwise_strip_params(nx, ny, c_lo, c_hi)
            wall_u = _cap_wall(
                skin_thickness=skin_thickness,
                x_start=x_lo,
                x_end=x_hi,
                cap_height=cap_h_upper,
                surface="suction_surface",
                fw=strip_c_fw,
                camber_xy=_camber_c,
            )
            wall_l = _cap_wall(
                skin_thickness=skin_thickness,
                x_start=x_lo,
                x_end=x_hi,
                cap_height=cap_h_lower,
                surface="pressure_surface",
                fw=strip_c_fw,
                camber_xy=_camber_c,
            )
            if wall_u is not None:
                recipes[f"spar_cap_suction_surface_{i}"] = wall_u
                cu = _project_cap(
                    wall_u,
                    airfoil_sdf,
                    airfoil_expr_node=airfoil_expr_node,
                    camber_xy=_camber_c,
                    strip_clip=strip_c,
                )
                caps_u.append(cu)
                unrotated[f"spar_cap_suction_surface_{i}"] = cu
            if wall_l is not None:
                recipes[f"spar_cap_pressure_surface_{i}"] = wall_l
                cl = _project_cap(
                    wall_l,
                    airfoil_sdf,
                    airfoil_expr_node=airfoil_expr_node,
                    camber_xy=_camber_c,
                    strip_clip=strip_c,
                )
                caps_l.append(cl)
                unrotated[f"spar_cap_pressure_surface_{i}"] = cl
        # Continuous labels point at the first surviving pad's recipe for
        # midline emit; discrete labels remain in recipes for constructive unions.
        for side, cont in (
            ("suction_surface", ComponentLabel.SPAR_CAP_SUCTION_SURFACE),
            ("pressure_surface", ComponentLabel.SPAR_CAP_PRESSURE_SURFACE),
        ):
            for i in range(N):
                key = f"spar_cap_{side}_{i}"
                if key in recipes:
                    recipes[_label_str(cont)] = recipes[key]
                    break
        if caps_u:
            cap_upper_expr = union_expr([c.expr for c in caps_u])
            if web_union is not None:
                cap_upper_expr = Subtract(cap_upper_expr, web_union)
            cap_upper = _wrap_expr_cap(caps_u[0], cap_upper_expr)
        if caps_l:
            cap_lower_expr = union_expr([c.expr for c in caps_l])
            if web_union is not None:
                cap_lower_expr = Subtract(cap_lower_expr, web_union)
            cap_lower = _wrap_expr_cap(caps_l[0], cap_lower_expr)

    return cap_upper, cap_lower, unrotated, recipes


__all__ = ["build_spar_caps"]
