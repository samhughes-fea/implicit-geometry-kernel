"""Section member builders: skin, shear webs, spar caps, sandwich cores, edge inserts.

Each builder turns a :class:`~.plan.SectionPlan` into ``structure`` contributions, in
the chord frame, and attaches the section component it made as the contribution's
``source`` (a :class:`~.plan.SectionMember`). Registration order is build order:
caps read the webs, and cores read the webs and caps.

How each member maps onto ``structure``:

* suction-side skin: the base wall; pressure-side skin: a wall partition on its side
  of the camber line (it sets the wall thickness there);
* spar caps: stacked wall layers beneath the skin, on their footprint;
* shear webs: branches whose mid-surface is the web axis;
* sandwich cores and edge inserts: regions.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, ClassVar, Protocol, runtime_checkable

import numpy as np

from implicit_geometry_kernel.contract.errors import InfeasibleSectionGeometry
from implicit_geometry_kernel.profiling import timed_stage
from implicit_geometry_kernel.sdf.nodes import HalfPlane
from implicit_geometry_kernel.sdf.nodes.tagged import ExprSDF
from implicit_geometry_kernel.section.assembly.build_spar_caps import build_spar_caps
from implicit_geometry_kernel.section.assembly.clip_algebra import (
    _parallel_web_half_gt,
    _parallel_web_half_lt,
    _parallel_web_strip_clip,
)
from implicit_geometry_kernel.section.assembly.labels import ComponentLabel
from implicit_geometry_kernel.section.assembly.plan import SectionMember, SectionPlan, b_frame_web_search
from implicit_geometry_kernel.section.assembly.web_layout import flapwise_web_sframe_inputs
from implicit_geometry_kernel.section.walls.walls import CoreWall, InsertWall, SkinWall, WebWall
from implicit_geometry_kernel.structure import (
    BranchContribution,
    ConstantLaw,
    Contribution,
    MemberBuilder,
    RegionContribution,
    WallContribution,
)


def _label_str(label) -> str:
    return str(getattr(label, "value", label))


def members_of(contributions: Sequence[Contribution], kind: str) -> list[SectionMember]:
    """The section members of the contributions whose feature kind is ``kind``."""
    return [c.source for c in contributions if c.spec.kind == kind and isinstance(c.source, SectionMember)]


@runtime_checkable
class WallType(Protocol):
    """One registered wall kind: labels, solids, midlines, and stage hooks."""

    kind: str
    has_midline: bool
    has_polygon: bool
    is_wall_solid: bool
    is_core: bool
    is_insert: bool
    cell_face_kind: str | None
    label_prefixes: ClassVar[tuple[str, ...]]
    exact_labels: ClassVar[tuple[str, ...]]

    def place(self, spec: Any) -> list[Any]: ...
    def field(self, inst: Any) -> Any: ...
    def polygon(self, inst: Any, boundary: Any) -> Any: ...
    def midline(self, inst: Any, boundary: Any) -> Any | None: ...
    def junctions(self, inst: Any, others: Any) -> list[Any]: ...


class WallTypeMixin:
    """Shared label matching and stage flags for :class:`SECTION_BUILDERS`."""

    kind: str = ""
    has_midline: bool = False
    has_polygon: bool = True
    is_wall_solid: bool = True
    is_core: bool = False
    is_insert: bool = False
    cell_face_kind: str | None = None
    label_prefixes: ClassVar[tuple[str, ...]] = ()
    exact_labels: ClassVar[tuple[str, ...]] = ()

    @classmethod
    def matches_label(cls, label: Any) -> bool:
        key = str(getattr(label, "value", label))
        if key in cls.exact_labels:
            return True
        return any(key.startswith(prefix) for prefix in cls.label_prefixes)

    @classmethod
    def format_label(cls, index: int) -> str:
        if cls.label_prefixes:
            prefix = cls.label_prefixes[0]
            return f"{prefix}{index}" if prefix.endswith("_") else f"{prefix}_{index}"
        return f"{cls.kind}_{index}"

    def place(self, spec: Any) -> list[Any]:
        return self.build(spec, ())  # type: ignore[attr-defined]

    def field(self, inst: Any) -> Any:
        return getattr(inst, "expr", None)

    def polygon(self, inst: Any, boundary: Any) -> Any:
        if not self.has_polygon:
            return None
        if hasattr(inst, "to_boundary_polygon"):
            return inst.to_boundary_polygon(boundary)
        return None

    def midline(self, inst: Any, boundary: Any) -> Any | None:
        if not self.has_midline:
            return None
        if hasattr(inst, "to_midline"):
            return inst.to_midline(boundary)
        return None

    def junctions(self, inst: Any, others: Any) -> list[Any]:
        return []


def wall_type_for_label(label: Any) -> type[WallTypeMixin] | None:
    """The registered builder whose labels match ``label``, or ``None``."""
    for builder in SECTION_BUILDERS:
        if builder.matches_label(label):
            return builder
    return None


def wall_cell_face(label: Any) -> str | None:
    wt = wall_type_for_label(label)
    return None if wt is None else wt.cell_face_kind


class SkinBuilder(WallTypeMixin, MemberBuilder):
    """Suction and pressure skins (camber-upper and camber-lower half annuli)."""

    kind = "skin"
    has_midline = True
    cell_face_kind = "skin"
    exact_labels = ("outer_skin",)
    label_prefixes = ("skin_suction_surface", "skin_pressure_surface")

    def build(self, ctx: SectionPlan, prior):  # type: ignore
        out: list[Contribution] = []
        with timed_stage("assemble_skin", ctx.stage_sink):
            for label, side, t in (
                (ComponentLabel.SKIN_SUCTION_SURFACE, "suction_surface", ctx.t_ss),
                (ComponentLabel.SKIN_PRESSURE_SURFACE, "pressure_surface", ctx.t_ps),
            ):
                wall = SkinWall(thickness=float(t), skin_side=side)
                comp = wall.to_sdf_component(ctx.airfoil_sdf, airfoil_expr_node=ctx.airfoil_expr)
                base = side == "suction_surface"
                out.append(WallContribution(
                    ctx.table.add(_label_str(label), "skin"),
                    ConstantLaw(t),
                    footprint=None if base else ExprSDF(comp.footprint_expr),
                    source=SectionMember(label, comp, {_label_str(label): wall}),
                ))
        return out


class WebBuilder(WallTypeMixin, MemberBuilder):
    """One shear web per anchored station."""

    kind = "web"
    has_midline = True
    cell_face_kind = "web"
    label_prefixes = ("shear_web_", "web_")

    def build(self, ctx: SectionPlan, prior):  # type: ignore
        out: list[Contribution] = []
        with timed_stage("assemble_webs", ctx.stage_sink):
            phi_inner_b, y_search_b = b_frame_web_search(ctx)
            for i, (x, t, align, (y_top, y_bot)) in enumerate(
                zip(ctx.xs, ctx.web_thicknesses, ctx.web_alignments, ctx.anchors)
            ):
                if ctx.explicit_endpoints is not None:
                    xt, yt, xb, yb = ctx.explicit_endpoints[i]
                elif align == "flapwise" and phi_inner_b is not None:
                    xt, yt, xb, yb = flapwise_web_sframe_inputs(
                        x, y_top, y_bot, ctx.twist, phi_inner_b, y_search_b
                    )
                else:
                    xt, yt, xb, yb = x, y_top, x, y_bot
                wall = WebWall(
                    skin_thickness=ctx.skin_thickness, x_top=float(xt), y_top=float(yt),
                    x_bot=float(xb), y_bot=float(yb), thickness=float(t), alignment=str(align),
                    twist_angle=ctx.twist,
                )
                comp = wall.to_sdf_component(ctx.airfoil_sdf, airfoil_expr_node=ctx.airfoil_expr)
                label = self.format_label(i)
                out.append(BranchContribution(
                    ctx.table.add(label, "web"),
                    mid=ExprSDF(_axis_line(xt, yt, xb, yb)),
                    region=ExprSDF(comp.expr),
                    thickness=ConstantLaw(t),
                    source=SectionMember(label, comp, {label: wall}),
                ))
        return out


def _axis_line(xt, yt, xb, yb) -> HalfPlane:
    """Signed distance to the infinite line through the web's two anchors."""
    dx, dy = float(xt) - float(xb), float(yt) - float(yb)
    length = float(np.hypot(dx, dy))
    ux, uy = (dx / length, dy / length) if length > 1e-12 else (0.0, 1.0)
    nx, ny = uy, -ux
    return HalfPlane(nx=nx, ny=ny, d=-(nx * float(xb) + ny * float(yb)))


class CapBuilder(WallTypeMixin, MemberBuilder):
    """Suction and pressure spar caps under the skin, placed by the cap rule."""

    kind = "cap"
    has_midline = True
    cell_face_kind = "cap"
    label_prefixes = ("spar_cap",)

    def build(self, ctx: SectionPlan, prior):  # type: ignore
        if ctx.cap_rule.kind == "none":
            return []
        webs = [m.component for m in members_of(prior, "web")]
        nx, ny, c_station = ctx.flapwise_stations()
        with timed_stage("assemble_caps", ctx.stage_sink):
            upper, lower, pads, recipes = build_spar_caps(
                ctx.airfoil_sdf, ctx.skin_thickness, ctx.xs, ctx.anchors,
                ctx.cap_h_upper, ctx.cap_h_lower, webs,
                cap_placement=ctx.cap_rule, airfoil_expr_node=ctx.airfoil_expr,
                all_flapwise=ctx.all_flapwise, c_station=c_station, nx=nx, ny=ny,
                N=ctx.n_webs, pitch_fr=ctx.pitch_fr, fhw=ctx.fhw, dhw=ctx.dhw,
                y_search=ctx.y_search, phi_inner_ref=ctx.phi_inner_ref,
                spar_cap_ply_pitch_m=ctx.spar_cap_ply_pitch_m,
            )
        # Run-out is decided before assembly (layout_without_spar_caps); a missing
        # side here would build a section with a hole where a cap belongs.
        if upper is None or lower is None:
            raise InfeasibleSectionGeometry(
                "Cap placement "
                f"{ctx.cap_rule.kind!r} expected spar caps but build_spar_caps "
                f"returned suction={upper is not None} "
                f"pressure={lower is not None}. Use "
                "layout_without_spar_caps when the cavity cannot hold "
                "minimum-gauge caps (run-out when a cap cannot be built at gauge); do not assemble None caps."
            )
        out: list[Contribution] = []
        for label, comp, height, first in (
            (ComponentLabel.SPAR_CAP_SUCTION_SURFACE, upper, ctx.cap_h_upper, True),
            (ComponentLabel.SPAR_CAP_PRESSURE_SURFACE, lower, ctx.cap_h_lower, False),
        ):
            out.append(WallContribution(
                ctx.table.add(_label_str(label), "cap"),
                ConstantLaw(height),
                footprint=ExprSDF(comp.footprint_expr),
                stacked=True,
                # The recipes and pads ride on the first cap so they are added once.
                source=SectionMember(label, comp, dict(recipes) if first else {}, dict(pads) if first else {}),
            ))
        return out


class CoreBuilder(WallTypeMixin, MemberBuilder):
    """One sandwich core per bay between webs, excluding the webs and caps."""

    kind = "core"
    has_midline = False
    is_wall_solid = False
    is_core = True
    cell_face_kind = None
    exact_labels = ("sandwich_core",)
    label_prefixes = ("sandwich_core_",)

    def build(self, ctx: SectionPlan, prior):  # type: ignore
        n = ctx.n_webs
        if not ctx.core_enabled:
            return []
        walls = members_of(prior, "web") if ctx.cap_rule.kind == "none" else (
            members_of(prior, "cap") + members_of(prior, "web"))
        exclusion = [m.component.expr for m in walls]
        nx, ny, c_station = ctx.flapwise_stations()
        bays = [(None, None)] if n == 0 else [
            (None if k == 0 else k - 1, None if k == n else k) for k in range(n + 1)
        ]
        out: list[Contribution] = []
        with timed_stage("assemble_cores", ctx.stage_sink):
            for k, (fore, aft) in enumerate(bays):
                if n > 0 and c_station is not None:
                    wall, sclip = _flapwise_core(ctx, nx, ny, c_station, fore, aft)
                    comp = wall.to_sdf_component(ctx.airfoil_sdf, airfoil_expr_node=ctx.airfoil_expr,
                                                 exclusion_exprs=exclusion, strip_clip=sclip)
                else:
                    wall = CoreWall(
                        skin_thickness=ctx.skin_thickness, bay_mode="chordwise",
                        x_start=None if fore is None else ctx.xs[fore],
                        x_end=None if aft is None else ctx.xs[aft],
                    )
                    comp = wall.to_sdf_component(ctx.airfoil_sdf, airfoil_expr_node=ctx.airfoil_expr,
                                                 exclusion_exprs=exclusion)
                label = self.format_label(k)
                out.append(RegionContribution(ctx.table.add(label, "core"), ExprSDF(comp.expr),
                                              source=SectionMember(label, comp, {label: wall})))
        return out


def _flapwise_core(ctx, nx, ny, c_station, fore, aft):
    if fore is None:
        sclip, mode, c_lo, c_hi = _parallel_web_half_lt(nx, ny, c_station[0]), "half_lt", float(c_station[0]), None
    elif aft is None:
        sclip, mode, c_lo, c_hi = _parallel_web_half_gt(nx, ny, c_station[-1]), "half_gt", float(c_station[-1]), None
    else:
        a, b = c_station[fore], c_station[aft]
        sclip, mode, c_lo, c_hi = _parallel_web_strip_clip(nx, ny, a, b), "strip", float(min(a, b)), float(max(a, b))
    wall = CoreWall(
        skin_thickness=ctx.skin_thickness, bay_mode="flapwise", flapwise_nx=float(nx),
        flapwise_ny=float(ny), flapwise_c_lo=c_lo, flapwise_c_hi=c_hi,
        flapwise_clip_mode=mode,  # type: ignore[arg-type]
    )
    return wall, sclip


class InsertBuilder(WallTypeMixin, MemberBuilder):
    """Trailing-edge and leading-edge inserts, when requested."""

    kind = "insert"
    has_midline = False
    is_insert = True
    cell_face_kind = None
    exact_labels = ("leading_edge_insert", "trailing_edge_insert")

    def build(self, ctx: SectionPlan, prior):  # type: ignore
        out: list[Contribution] = []
        if ctx.trailing_edge_insert_x_m is not None:
            wall = InsertWall(kind="trailing_edge", skin_thickness=ctx.skin_thickness,
                              x_limit=float(ctx.trailing_edge_insert_x_m))
            out.append(self._insert(ctx, ComponentLabel.TRAILING_EDGE_INSERT, wall))
        if ctx.leading_edge_insert_x_m is not None:
            try:
                le_x = float(ctx.airfoil_sdf.leading_edge[0])
            except AttributeError:
                le_x = 0.0
            wall = InsertWall(kind="leading_edge", skin_thickness=ctx.skin_thickness,
                              x_limit=float(ctx.leading_edge_insert_x_m), le_x=le_x, le_y=0.0,
                              radius=ctx.le_radius)
            out.append(self._insert(ctx, ComponentLabel.LEADING_EDGE_INSERT, wall))
        return out

    @staticmethod
    def _insert(ctx, label, wall) -> RegionContribution:
        comp = wall.to_sdf_component(ctx.airfoil_sdf, airfoil_expr_node=ctx.airfoil_expr)
        return RegionContribution(ctx.table.add(_label_str(label), "insert"), ExprSDF(comp.expr),
                                  source=SectionMember(label, comp, {_label_str(label): wall}))


SECTION_BUILDERS = (SkinBuilder, WebBuilder, CapBuilder, CoreBuilder, InsertBuilder)
"""Default build order for a section."""

__all__ = [
    "CapBuilder",
    "CoreBuilder",
    "InsertBuilder",
    "SECTION_BUILDERS",
    "SkinBuilder",
    "WallType",
    "WallTypeMixin",
    "WebBuilder",
    "members_of",
    "wall_cell_face",
    "wall_type_for_label",
]
