"""per_web_pads spar-cap branch (peeled from build_spar_caps)."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from implicit_geometry_kernel.sdf.nodes import Expr, Subtract, union_expr
from implicit_geometry_kernel.section.assembly.cap_feasibility import cap_feasibility_clamp
from implicit_geometry_kernel.section.assembly.clip_algebra import (
    _parallel_web_strip_clip,
    flapwise_strip_params,
    web_station_projection,
)
from implicit_geometry_kernel.section.assembly.labels import ComponentLabel
from implicit_geometry_kernel.section.walls.spar_cap import SparCap
from implicit_geometry_kernel.section.walls.walls import CapWall

_ONE_PLY_REL_TOL = 1e-9


def snap_cap_height_at_most(height_m: float, pitch_m: float | None) -> float:
    """Greatest on-grid height that does not exceed `height_m`.

    Returns `height_m` unchanged when `pitch_m` is absent/invalid.
    When not even one ply fits (`height_m < pitch`), returns ``0.0`` — geometric
    run-out, matching ``PlyGrid.largest_at_most``. Leaving a sub-ply continuous
    strip made the ply-count ceil (1 ply) disagree with geometry and broke Catalog FI.
    """
    if pitch_m is None:
        return float(height_m)
    p = float(pitch_m)
    c = float(height_m)
    if not (p > 0.0) or not math.isfinite(p) or not math.isfinite(c):
        return float(height_m)
    if c < p:
        return 0.0
    return math.floor(c / p + _ONE_PLY_REL_TOL) * p


def snap_cap_heights(
    h_upper: float,
    h_lower: float,
    pitch_m: float | None,
) -> tuple[float, float]:
    """Snap a suction/pressure pair after feasibility clamp."""
    return (
        snap_cap_height_at_most(h_upper, pitch_m),
        snap_cap_height_at_most(h_lower, pitch_m),
    )


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
    cap = wall.to_sdf_component(airfoil_sdf, camber_xy=camber_xy, airfoil_expr_node=airfoil_expr_node)
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


def fill_per_web_pads(
    *,
    airfoil_sdf,
    skin_thickness,
    xs,
    anchors,
    cap_h_upper,
    cap_h_lower,
    N: int,
    dhw: float,
    y_search,
    phi_inner_ref,
    spar_cap_ply_pitch_m,
    all_flapwise: bool,
    c_station,
    nx,
    ny,
    recipes: dict,
    unrotated: dict,
    airfoil_expr_node,
    web_union,
    _cap_wall,
    _project_cap,
    _wrap_expr_cap,
    _label_str,
) -> tuple[Any | None, Any | None]:
    if N < 1:
        raise ValueError(f"per_web_pads cap placement requires at least one web (got N={N}).")
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
        cap_h_upper = min(float(c[2]) for c in _clamped_c)
        cap_h_lower = min(float(c[3]) for c in _clamped_c)
    cap_h_upper, cap_h_lower = snap_cap_heights(cap_h_upper, cap_h_lower, spar_cap_ply_pitch_m)
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
    cap_upper = cap_lower = None
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
            strip_c_fw = flapwise_strip_params(nx, ny, c_lo, c_hi)
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
        su_key = f"spar_cap_suction_surface_{i}"
        ps_key = f"spar_cap_pressure_surface_{i}"
        if wall_u is not None:
            recipes[su_key] = wall_u
            cu = _project_cap(
                wall_u,
                airfoil_sdf,
                airfoil_expr_node=airfoil_expr_node,
                camber_xy=_camber_c,
                strip_clip=strip_c,
            )
            caps_u.append(cu)
            unrotated[su_key] = cu
        if wall_l is not None:
            recipes[ps_key] = wall_l
            cl = _project_cap(
                wall_l,
                airfoil_sdf,
                airfoil_expr_node=airfoil_expr_node,
                camber_xy=_camber_c,
                strip_clip=strip_c,
            )
            caps_l.append(cl)
            unrotated[ps_key] = cl
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
        footprint = union_expr([c.footprint_expr for c in caps_u])
        if web_union is not None:
            cap_upper_expr = Subtract(cap_upper_expr, web_union)
            footprint = Subtract(footprint, web_union)
        cap_upper = _wrap_expr_cap(caps_u[0], cap_upper_expr)
        cap_upper.footprint_expr = footprint
    if caps_l:
        cap_lower_expr = union_expr([c.expr for c in caps_l])
        footprint = union_expr([c.footprint_expr for c in caps_l])
        if web_union is not None:
            cap_lower_expr = Subtract(cap_lower_expr, web_union)
            footprint = Subtract(footprint, web_union)
        cap_lower = _wrap_expr_cap(caps_l[0], cap_lower_expr)
        cap_lower.footprint_expr = footprint
    return cap_upper, cap_lower
