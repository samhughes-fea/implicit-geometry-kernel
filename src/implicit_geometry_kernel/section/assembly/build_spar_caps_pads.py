"""per_web_pads spar-cap branch (peeled from build_spar_caps)."""

from __future__ import annotations

from typing import Any

from implicit_geometry_kernel.sdf.nodes import Subtract, union_expr
from implicit_geometry_kernel.section.assembly.cap_feasibility import cap_feasibility_clamp
from implicit_geometry_kernel.section.assembly.cap_ply_snap import snap_cap_heights
from implicit_geometry_kernel.section.assembly.clip_algebra import (
    _parallel_web_strip_clip,
    flapwise_strip_params,
    web_station_projection,
)
from implicit_geometry_kernel.section.assembly.labels import ComponentLabel


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
        raise ValueError(
            "per_web_pads cap placement requires at least one web "
            f"(got N={N})."
        )
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
    cap_h_upper, cap_h_lower = snap_cap_heights(
        cap_h_upper, cap_h_lower, spar_cap_ply_pitch_m
    )
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
        su_key = "spar_cap_suction_surface" if i == 0 else f"spar_cap_suction_surface_{i}"
        ps_key = "spar_cap_pressure_surface" if i == 0 else f"spar_cap_pressure_surface_{i}"
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
            key = f"spar_cap_{side}" if i == 0 else f"spar_cap_{side}_{i}"
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
