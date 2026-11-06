"""Sections are assembled through ``structure``: check the member model they carry."""

import numpy as np
import pytest

from implicit_geometry_kernel.sdf.airfoil import AirfoilSDF
from implicit_geometry_kernel.sdf.ir import Circle
from implicit_geometry_kernel.sdf.nodes import ExprSDF
from implicit_geometry_kernel.section.api.factory import build_section_view, resolve_system_type
from implicit_geometry_kernel.section.assembly.members import SECTION_BUILDERS
from implicit_geometry_kernel.section.assembly.multicell import MultiCellSection, section_registry
from implicit_geometry_kernel.section.assembly.plan import SectionMember
from implicit_geometry_kernel.structure import MemberBuilder, RegionContribution

SKIN, CAP, WEB = 0.004, 0.012, 0.005


def host(key, twist=0.0):
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    view = build_section_view(airfoil, resolve_system_type(key), twist_angle_rad=twist,
                              skin_thickness=SKIN, spar_cap_thickness=CAP, shear_web_thickness=WEB)
    section = view
    while getattr(section, "underlying_section", section) is not section:
        section = section.underlying_section
    return section


def owner(structure, x, y):
    s = structure.material.evaluate(np.array([[x, y, 0.0]]))
    return structure.table[int(s.meta[0])].name if s.phi[0] < 0 else None


@pytest.mark.parametrize("key", ["0A", "0B", "1C-CN", "2D-F", "3B-CN"])
def test_every_component_is_a_member(key):
    section = host(key)
    names = {spec.name for spec in section.structure.table}
    labels = {str(getattr(lbl, "value", lbl)) for lbl in section.labels}
    assert labels == names


def test_box_section_members_and_fields():
    section = host("2D-CN")
    s = section.structure
    kinds = {spec.name: spec.kind for spec in s.table}
    assert kinds["skin_suction_surface"] == "skin" and kinds["skin_pressure_surface"] == "skin"
    assert kinds["spar_cap_suction_surface"] == "cap" and kinds["shear_web_0"] == "web"
    assert [b.spec.name for b in s.branches] == ["shear_web_0", "shear_web_1"]
    assert all(w.stacked for w in s.walls if w.spec.kind == "cap")

    # Find the upper surface above the box mid-point, then probe down through the wall.
    x = 0.5 * (s.branches[0].contribution.source.recipes["shear_web_0"].x_top
               + s.branches[1].contribution.source.recipes["shear_web_1"].x_top)
    ys = np.linspace(0.0, 0.2, 20001)
    phi = s.oml.evaluate(np.column_stack([np.full_like(ys, x), ys, np.zeros_like(ys)])).phi
    y_top = ys[np.argmin(np.abs(phi))]
    assert s.thickness(np.array([[x, y_top, 0.0]]))[0] == pytest.approx(SKIN + CAP)
    assert owner(s, x, y_top - SKIN / 2) == "skin_suction_surface"
    assert owner(s, x, y_top - SKIN - CAP / 2) == "spar_cap_suction_surface"
    assert owner(s, x, y_top - SKIN - CAP - 0.005) == "sandwich_core_1"


def test_webs_are_branches_on_their_axis():
    s = host("2D-CN").structure
    for b in s.branches:
        wall = next(iter(b.contribution.source.recipes.values()))
        on_axis = np.array([[wall.x_top, 0.5 * (wall.y_top + wall.y_bot), 0.0]])
        assert b.contribution.mid.evaluate(on_axis).phi[0] == pytest.approx(0.0, abs=1e-12)
        assert b.mask.evaluate(on_axis).phi[0] < 0.0  # inside the wall reference surface
        assert owner(s, float(on_axis[0, 0]), float(on_axis[0, 1])) == b.spec.name


def test_a_new_member_type_is_one_more_builder():
    class Ballast(MemberBuilder):
        def build(self, ctx, prior):
            comp = ExprSDF(Circle(cx=0.6, cy=0.0, r=0.01))
            spec = ctx.table.add("ballast", "bond")
            return [RegionContribution(spec, comp, source=None)]

    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    registry = section_registry().register(Ballast())
    section = MultiCellSection(airfoil, [0.25, 0.5], registry=registry, skin_thickness=SKIN)
    assert owner(section.structure, 0.6, 0.0) == "ballast"
    assert "ballast" not in {str(getattr(lbl, "value", lbl)) for lbl in section.labels}


def test_default_builders_and_sources():
    assert [b.__name__ for b in SECTION_BUILDERS] == [
        "SkinBuilder", "WebBuilder", "CapBuilder", "CoreBuilder", "InsertBuilder"]
    section = host("1C-CN", twist=0.3)
    sources = [c.source for c in (*section.structure.walls, *section.structure.regions)]
    assert all(isinstance(m, SectionMember) for m in sources)
