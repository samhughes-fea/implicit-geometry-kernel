"""Everything decided before a section's members are built.

:func:`plan_section` validates the inputs of :class:`MultiCellSection`, resolves the
cap placement rule, clamps the gauges and anchors the webs on the inner skin. The
result, :class:`SectionPlan`, is the :class:`~implicit_geometry_kernel.structure.BuildContext`
that the member builders (skin, webs, caps, cores, inserts) read.

The section's structure is built in the chord frame; twist is applied afterwards.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from implicit_geometry_kernel.contract.errors import InfeasibleSectionGeometry
from implicit_geometry_kernel.profiling import timed_stage
from implicit_geometry_kernel.sdf.nodes import Offset
from implicit_geometry_kernel.sdf.nodes.tagged import ExprSDF
from implicit_geometry_kernel.section.assembly.clip_algebra import airfoil_expr, flapwise_web_normal, web_station_projection
from implicit_geometry_kernel.section.assembly.web_layout import resolve_web_anchors
from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
from implicit_geometry_kernel.section.layout.structural import parse_fixed_cap_anchor, parse_structural_family
from implicit_geometry_kernel.section.layout.thickness_limits import clamp_skin_thickness_m, clamp_spar_thickness_m, clamp_web_thickness_m
from implicit_geometry_kernel.section.walls.base import expr_callable
from implicit_geometry_kernel.structure import BuildContext


@dataclass(frozen=True)
class SectionMember:
    """What a section builder attaches to each contribution as its ``source``.

    ``key`` is the component label, ``component`` the chord-frame component,
    ``recipes`` the wall recipes it adds (label to recipe) and ``extra`` any further
    chord-frame components (the discrete cap pads).
    """

    key: Any
    component: Any
    recipes: dict[str, Any] = field(default_factory=dict)
    extra: dict[Any, Any] = field(default_factory=dict)


@dataclass(frozen=True, kw_only=True)
class SectionPlan(BuildContext):
    """Resolved inputs for one section (chord frame, metres, radians)."""

    airfoil_sdf: Any
    airfoil_expr: Any
    twist: float
    cap_rule: CapPlacementRule
    structural_family: str
    t_ss: float
    t_ps: float
    skin_thickness: float
    phi_inner_ref: Any
    xs: list[float]
    web_indices: tuple[int, ...]
    anchors: list[tuple[float, float]]
    explicit_endpoints: list[tuple[float, float, float, float]] | None
    web_thicknesses: list[float]
    web_alignments: list[str]
    cap_h_upper: float
    cap_h_lower: float
    pitch_fr: float
    fhw: float
    dhw: float
    y_search: Any
    core_enabled: bool
    trailing_edge_insert_x_m: float | None
    leading_edge_insert_x_m: float | None
    le_radius: float | None
    spar_cap_ply_pitch_m: float | None
    stage_sink: dict[str, float] | None

    @property
    def n_webs(self) -> int:
        return len(self.xs)

    @property
    def all_flapwise(self) -> bool:
        return self.n_webs > 0 and all(a == "flapwise" for a in self.web_alignments)

    def flapwise_stations(self):
        """``(nx, ny, c_station)`` for all-flapwise webs, else ``(None, None, None)``."""
        if not self.all_flapwise:
            return None, None, None
        nx, ny = flapwise_web_normal(float(self.twist))
        c = [web_station_projection(nx, ny, self.xs[i], 0.5 * (self.anchors[i][0] + self.anchors[i][1]))
             for i in range(self.n_webs)]
        return nx, ny, c


def _cap_rule(structural_family, cap_placement, fixed_cap_anchor) -> tuple[str, CapPlacementRule]:
    fam = parse_structural_family(structural_family)
    if cap_placement is not None:
        return fam, cap_placement
    pos = parse_fixed_cap_anchor(fixed_cap_anchor) if fam == "B" else None
    return fam, CapPlacementRule.from_structural_family(fam, spar_cap_positioning=pos)


def _per_web(value, n: int, name: str, cast) -> list:
    if isinstance(value, (int, float, str)):
        return [cast(value)] * n
    out = [cast(v) for v in value]
    if len(out) != n:
        raise ValueError(f"{name} must be scalar or length {n} (got {len(out)}).")
    return out


def _anchor_webs(airfoil_sdf, skin_thickness, xs_requested, *, explicit_endpoints, web_y_coords,
                 y_search, phi_inner_ref):
    if explicit_endpoints is not None:
        return ([0.5 * (e[0] + e[2]) for e in explicit_endpoints],
                [(e[1], e[3]) for e in explicit_endpoints])
    if web_y_coords is not None:
        # Frozen anchors from resolve_station_layout: never re-probe (that could
        # shorten the web set).
        xs = list(xs_requested)
        anchors = [(float(yt), float(yb)) for yt, yb in web_y_coords]
        if len(anchors) != len(xs):
            raise ValueError(f"web_y_coords must have {len(xs)} entries.")
        return xs, anchors
    # Fail closed: every requested web must anchor. Stepping the layout down
    # (2D-F -> 1C-F) is resolve_station_layout's job, before assembly.
    resolved, dropped = resolve_web_anchors(
        airfoil_sdf, skin_thickness, xs_requested, y_search=y_search,
        phi_inner=phi_inner_ref, return_dropped=True,
    )
    if xs_requested and not resolved:
        raise InfeasibleSectionGeometry(
            "No feasible web anchors: local thickness is smaller than "
            f"inner-skin offset (skin_thickness={skin_thickness:.6g}) "
            f"at requested x={xs_requested}."
        )
    if dropped:
        drop_desc = ", ".join(f"{float(d.x_m):.6g} ({d.reason})" for d in dropped)
        raise InfeasibleSectionGeometry(
            "Requested webs cannot all anchor at this skin_thickness="
            f"{skin_thickness:.6g}: dropped [{drop_desc}]. "
            "Call resolve_station_layout (or pass frozen web_y_coords) "
            "so the layout steps down before assemble; assemble does "
            "not silently shrink the web set."
        )
    return [a[0] for a in resolved], [(a[1], a[2]) for a in resolved]


def plan_section(
    airfoil_sdf, web_x_positions, web_thickness, web_alignment, cap_height, skin_thickness,
    twist_angle, web_y_coords, trailing_edge_insert_x_m, leading_edge_insert_x_m, le_radius,
    core_enabled, y_search, structural_family, fixed_cap_anchor, pitch_fraction_of_chord_from_le,
    fixed_cap_chord_half_width, discrete_cap_chord_half_width, *, cap_placement, web_endpoints,
    stage_sink, spar_cap_ply_pitch_m, web_indices=None,
) -> SectionPlan:
    """Validate and resolve the inputs of :class:`MultiCellSection` (same arguments)."""
    explicit = None
    if web_endpoints is not None:
        explicit = [(float(xt), float(yt), float(xb), float(yb)) for xt, yt, xb, yb in web_endpoints]
        xs_requested = [0.5 * (e[0] + e[2]) for e in explicit]
        idx_requested = list(range(len(xs_requested))) if web_indices is None else list(web_indices)
    else:
        xs_requested = [float(x) for x in web_x_positions]
        idx_requested = list(range(len(xs_requested))) if web_indices is None else [int(i) for i in web_indices]
        if len(idx_requested) != len(xs_requested):
            idx_requested = list(range(len(xs_requested)))
        order = sorted(range(len(xs_requested)), key=lambda i: xs_requested[i])
        xs_requested = [xs_requested[i] for i in order]
        idx_requested = [idx_requested[i] for i in order]
    family, rule = _cap_rule(structural_family, cap_placement, fixed_cap_anchor)
    if len(xs_requested) < 1 and rule.kind not in ("none", "fixed_band"):
        raise ValueError(
            "web_x_positions is empty, which is only supported for "
            f"cap placement 'none' or 'fixed_band' (got {rule.kind!r}). "
            "per_web_pads needs one pad per web; continuous_box spans "
            "outermost webs."
        )
    if isinstance(skin_thickness, (int, float)):
        t_ss = t_ps = clamp_skin_thickness_m(skin_thickness)
    else:
        t_ss = clamp_skin_thickness_m(float(skin_thickness[0]))
        t_ps = clamp_skin_thickness_m(float(skin_thickness[1]))
    t_skin = max(t_ss, t_ps)
    if rule.kind in ("fixed_band", "per_web_pads"):
        from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF

        if not isinstance(airfoil_sdf, AirfoilSDF):
            raise TypeError("fixed_band / per_web_pads require airfoil_sdf "
                            "to be AirfoilSDF (pitching / max-thickness anchors).")
    poly = airfoil_expr(airfoil_sdf)
    if poly is None:
        raise ValueError("MultiCellSection requires AirfoilSDF with vertices.")

    chord = float(getattr(airfoil_sdf, "chord", 1.0))
    fhw = (float(fixed_cap_chord_half_width) if fixed_cap_chord_half_width is not None
           else float(rule.half_width) if rule.half_width is not None else 0.05 * chord)
    if rule.kind == "fixed_band" and rule.half_width is None:
        rule = CapPlacementRule.fixed_band(anchor=rule.anchor or "max_thickness", half_width=fhw)
    dhw = float(discrete_cap_chord_half_width) if discrete_cap_chord_half_width is not None else 0.04 * chord

    phi_inner_ref = expr_callable(Offset(poly, amount=-t_skin))
    with timed_stage("assemble_anchors", stage_sink):
        xs, anchors = _anchor_webs(airfoil_sdf, t_skin, xs_requested, explicit_endpoints=explicit,
                                   web_y_coords=web_y_coords, y_search=y_search, phi_inner_ref=phi_inner_ref)
    n = len(xs)
    if rule.kind == "continuous_box" and n < 2:
        raise InfeasibleSectionGeometry(
            f"continuous_box cap placement requires at least 2 webs (got N={n}). "
            "There is no 1D catalogue layout; resolve_station_layout must step "
            "D-family to 1C (per_web_pads) before assemble."
        )
    web_t = [clamp_web_thickness_m(t) for t in _per_web(web_thickness, n, "web_thickness", float)]
    aligns = _per_web(web_alignment, n, "web_alignment", str)
    if isinstance(cap_height, (int, float)):
        cap_u = cap_l = float(cap_height)
    else:
        cap_u, cap_l = float(cap_height[0]), float(cap_height[1])

    return SectionPlan(
        oml=ExprSDF(Offset(poly, amount=0.0)),
        span=(0.0, 0.0),
        airfoil_sdf=airfoil_sdf,
        airfoil_expr=poly,
        twist=float(twist_angle),
        cap_rule=rule,
        structural_family=family,
        t_ss=float(t_ss),
        t_ps=float(t_ps),
        skin_thickness=float(t_skin),
        phi_inner_ref=phi_inner_ref,
        xs=xs,
        web_indices=tuple(idx_requested[: len(xs)]),
        anchors=anchors,
        explicit_endpoints=explicit,
        web_thicknesses=web_t,
        web_alignments=aligns,
        cap_h_upper=clamp_spar_thickness_m(cap_u),
        cap_h_lower=clamp_spar_thickness_m(cap_l),
        pitch_fr=float(pitch_fraction_of_chord_from_le),
        fhw=fhw,
        dhw=dhw,
        y_search=y_search,
        core_enabled=bool(core_enabled),
        trailing_edge_insert_x_m=trailing_edge_insert_x_m,
        leading_edge_insert_x_m=leading_edge_insert_x_m,
        le_radius=le_radius,
        spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
        stage_sink=stage_sink,
    )


def b_frame_web_search(plan: SectionPlan):
    """``(phi_inner_b, y_search_b)`` for flapwise webs on a twisted section, else ``(None, y_search)``."""
    from implicit_geometry_kernel.sdf.kernels.transforms import rotate_field

    flapwise = any(a == "flapwise" for a in plan.web_alignments)
    if not (flapwise and abs(plan.twist) > 1e-10):
        return None, plan.y_search
    phi_inner_b = rotate_field(plan.phi_inner_ref, plan.twist)
    y_search_b = plan.y_search
    verts = getattr(plan.airfoil_sdf, "vertices", None)
    if verts is not None:
        arr = np.asarray(verts, dtype=float)
        c, s = float(np.cos(plan.twist)), float(np.sin(plan.twist))
        yb = s * arr[:, 0] + c * arr[:, 1]
        y_search_b = (float(np.min(yb)) - 0.1, float(np.max(yb)) + 0.1)
    elif y_search_b is None:
        y_search_b = (-1.5, 1.5)
    return phi_inner_b, y_search_b


__all__ = ["SectionMember", "SectionPlan", "b_frame_web_search", "plan_section"]
