"""A dummy wall type registered only in tests flows through section stages."""

from __future__ import annotations

from implicit_geometry_kernel import Airfoil
from implicit_geometry_kernel.sdf.nodes import HalfPlane
from implicit_geometry_kernel.sdf.nodes.tagged import TagNode
from implicit_geometry_kernel.section import SectionOptions, SectionSpec, resolve_system_type
from implicit_geometry_kernel.section.assembly.members import (
    SECTION_BUILDERS,
    BondlineBuilder,
    CapBuilder,
    WallTypeMixin,
    wall_type_for_label,
)
from implicit_geometry_kernel.sdf.nodes import Circle
from implicit_geometry_kernel.section.walls.base import SectionComponent
from implicit_geometry_kernel.structure import RegionContribution
from implicit_geometry_kernel.section.assembly.plan import SectionMember
from implicit_geometry_kernel.section.spec.build import build_section_geometry
from implicit_geometry_kernel.section.walls.walls import StiffenerWall
from implicit_geometry_kernel.structure import BranchContribution, ConstantLaw, MemberBuilder
import numpy as np


class DummyStiffenerBuilder(WallTypeMixin, MemberBuilder):
    kind = "stiffener"
    has_midline = True
    on_shell_product = True
    is_wall_solid = True
    cell_face_kind = "stiffener"
    label_prefixes = ("stiffener_",)

    def build(self, ctx, prior):
        x = float(ctx.xs[0]) if ctx.xs else 0.35 * float(ctx.airfoil_sdf.chord or 1.0)
        y_top, y_bot = ctx.anchors[0] if ctx.anchors else (0.1, -0.1)
        t = float(ctx.web_thicknesses[0] if ctx.web_thicknesses else 0.004)
        wall = StiffenerWall(
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
                mid=TagNode(mid, 0),
                region=TagNode(comp.expr, 0),
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
        resolve_system_type("2D-F-LC-0R"),
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


def test_bondline_region_builder_is_registered_and_idle_by_default():
    assert wall_type_for_label("bondline_0") is BondlineBuilder
    assert BondlineBuilder in SECTION_BUILDERS
    spec = SectionSpec.from_kwargs(
        resolve_system_type("2D-F-LC-0R"),
        Airfoil.naca4("2412").to_sdf(chord=1.0),
        skin_thickness=0.004,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
    )
    result = build_section_geometry(
        spec,
        options=SectionOptions(outputs=("properties", "midlines", "cell_wall_loops")),
    )
    assert all("bondline" not in str(k) for k in (result.properties_map or ()))


def test_bondline_region_builder_emits_only_when_build_is_overridden(monkeypatch):
    def build(self, ctx, prior):
        label = self.format_label(0)
        comp = SectionComponent(Circle(cx=0.4, cy=0.0, r=0.002))
        return [
            RegionContribution(
                ctx.table.add(label, "bondline"),
                TagNode(comp.expr, 0),
                source=SectionMember(label, comp),
            )
        ]

    monkeypatch.setattr(BondlineBuilder, "build", build)
    spec = SectionSpec.from_kwargs(
        resolve_system_type("2D-F-LC-0R"),
        Airfoil.naca4("2412").to_sdf(chord=1.0),
        skin_thickness=0.004,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
    )
    result = build_section_geometry(
        spec,
        options=SectionOptions(outputs=("properties",)),
    )
    assert "bondline_0" in (result.properties_map or {})


def test_cap_builder_is_off_the_shell_product() -> None:
    assert CapBuilder.on_shell_product is False
    assert CapBuilder.polygon_by_surface is True
    assert DummyStiffenerBuilder.on_shell_product is True
    assert wall_type_for_label("stiffener_0") is None


def test_section_modules_do_not_hardcode_wall_prefixes() -> None:
    import ast
    from pathlib import Path

    banned = ("shear_web_", "spar_cap")
    root = Path(__file__).resolve().parents[2] / "src" / "implicit_geometry_kernel" / "section"
    hits: list[str] = []
    for path in sorted(root.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            name = getattr(func, "attr", None) or getattr(func, "id", None)
            if name != "startswith" or not node.args:
                continue
            arg = node.args[0]
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                if any(arg.value.startswith(p) for p in banned):
                    hits.append(f"{path.relative_to(root)}:{node.lineno}:{arg.value!r}")
    assert hits == []
