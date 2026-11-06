"""A dummy wall type registered only in tests flows through section stages."""

from __future__ import annotations

from implicit_geometry_kernel import Airfoil
from implicit_geometry_kernel.sdf.nodes import HalfPlane
from implicit_geometry_kernel.sdf.nodes.tagged import ExprSDF
from implicit_geometry_kernel.section import SectionOptions, SectionSpec, resolve_system_type
from implicit_geometry_kernel.section.assembly.members import (
    SECTION_BUILDERS,
    WallTypeMixin,
    wall_type_for_label,
)
from implicit_geometry_kernel.section.assembly.plan import SectionMember
from implicit_geometry_kernel.section.spec.build import build_section_geometry
from implicit_geometry_kernel.section.walls.walls import WebWall
from implicit_geometry_kernel.structure import BranchContribution, ConstantLaw, MemberBuilder
import numpy as np


class DummyStiffenerBuilder(WallTypeMixin, MemberBuilder):
    kind = "stiffener"
    has_midline = True
    is_wall_solid = True
    cell_face_kind = None
    label_prefixes = ("stiffener_",)

    def build(self, ctx, prior):
        x = float(ctx.xs[0]) if ctx.xs else 0.35 * float(ctx.airfoil_sdf.chord or 1.0)
        y_top, y_bot = ctx.anchors[0] if ctx.anchors else (0.1, -0.1)
        t = float(ctx.web_thicknesses[0] if ctx.web_thicknesses else 0.004)
        wall = WebWall(
            skin_thickness=ctx.skin_thickness,
            x_top=x,
            y_top=float(y_top),
            x_bot=x,
            y_bot=float(y_bot),
            thickness=t,
            alignment="chord_normal",
            twist_angle=ctx.twist,
        )
        comp = wall.to_sdf_component(ctx.airfoil_sdf, airfoil_expr_node=ctx.airfoil_expr)
        label = self.format_label(0)
        dx, dy = 0.0, float(y_top) - float(y_bot)
        length = float(np.hypot(dx, dy)) or 1.0
        nx, ny = dy / length, -dx / length
        mid = HalfPlane(nx=nx, ny=ny, d=-(nx * x + ny * float(y_bot)))
        return [
            BranchContribution(
                ctx.table.add(label, "stiffener"),
                mid=ExprSDF(mid),
                region=ExprSDF(comp.expr),
                thickness=ConstantLaw(t),
                source=SectionMember(label, comp, {label: wall}),
            )
        ]


def test_dummy_stiffener_flows_through_properties_midlines_and_cells(monkeypatch):
    monkeypatch.setattr(
        "implicit_geometry_kernel.section.assembly.members.SECTION_BUILDERS",
        SECTION_BUILDERS + (DummyStiffenerBuilder,),
    )
    assert wall_type_for_label("stiffener_0") is DummyStiffenerBuilder

    spec = SectionSpec.from_kwargs(
        resolve_system_type("2D-F"),
        Airfoil.naca4("2412").to_sdf(chord=1.0),
        skin_thickness=0.004,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
    )
    result = build_section_geometry(
        spec,
        options=SectionOptions(outputs=("properties", "midlines", "cell_wall_loops")),
    )
    assert "stiffener_0" in result.properties_map
    assert any(s.label == "stiffener_0" for s in result.midline_strips)
    assert result.cell_wall_loops
