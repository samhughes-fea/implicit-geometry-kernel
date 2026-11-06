"""A box beam built from all four contribution kinds, checked point by point."""

from collections.abc import Sequence

import numpy as np
import pytest

from implicit_geometry_kernel.sdf.nodes.tagged import PolygonSDF
from implicit_geometry_kernel.sdf.nodes.span import ComplementNode, MaxNode
from implicit_geometry_kernel.structure import (
    BranchContribution,
    BuildContext,
    BuilderRegistry,
    ConstantLaw,
    Contribution,
    MemberBuilder,
    RegionContribution,
    StructureAssembler,
    StructureError,
    TableLaw,
    TransverseContribution,
    WallContribution,
    assemble,
)

W, H, L = 1.0, 0.2, 6.0  # box width, height, span
SKIN, CAP, WEB = 0.01, 0.03, 0.008


def rect(x0, x1, y0, y1):
    return PolygonSDF(np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]]))


class Skin(MemberBuilder):
    def build(self, ctx, prior):
        return [WallContribution(ctx.table.add("skin", "skin"), ConstantLaw(SKIN))]


class Caps(MemberBuilder):
    def build(self, ctx, prior):
        top = rect(-0.2, 0.2, 0.0, 1.0)  # footprint: upper wall between x = -0.2 and 0.2
        return [WallContribution(ctx.table.add("cap_upper", "cap"), ConstantLaw(CAP), footprint=top)]


class Web(MemberBuilder):
    def build(self, ctx, prior):
        mid = rect(-1e-9, 1e-9, -1.0, 1.0)
        region = rect(-WEB / 2, WEB / 2, -1.0, 1.0)
        return [BranchContribution(ctx.table.add("web", "web"), mid, region, ConstantLaw(WEB))]


class Stringer(MemberBuilder):
    """Stands 0.04 m off the lower wall at x = -0.35, from z = 1 to z = 4."""

    def build(self, ctx, prior):
        region = rect(-0.355, -0.345, -1.0, -H / 2 + 0.04)
        mid = rect(-0.35 - 1e-9, -0.35 + 1e-9, -1.0, 1.0)
        lower = ComplementNode(rect(-1.0, 1.0, 0.0, 1.0))  # negative below the chord line
        return [BranchContribution(ctx.table.add("stringer", "stringer"), mid, region, ConstantLaw(0.01),
                                   clip=lower, depth=ConstantLaw(0.04), free_end=True, extent=(1.0, 4.0))]


class Core(MemberBuilder):
    def build(self, ctx, prior):
        bay = MaxNode(rect(0.25, 0.45, -1.0, 1.0), ctx.oml)
        return [RegionContribution(ctx.table.add("core_aft", "core"), bay)]


class Rib(MemberBuilder):
    def build(self, ctx, prior):
        return [TransverseContribution(ctx.table.add("rib_3", "transverse"), 3.0, ctx.oml, 0.02)]


@pytest.fixture(scope="module")
def structure():
    ctx = BuildContext(oml=rect(-W / 2, W / 2, -H / 2, H / 2), span=(0.0, L))
    registry = BuilderRegistry()
    for b in (Skin(), Caps(), Web(), Stringer(), Core(), Rib()):
        registry.register(b)
    return StructureAssembler(registry).assemble(ctx)


def owner(structure, x, y, z):
    s = structure.material.evaluate(np.array([[x, y, z]]))
    return (structure.table[int(s.meta[0])].name if s.phi[0] < 0 else None)


def test_wall_thickness_is_raised_inside_the_cap_footprint(structure):
    pts = np.array([[0.0, H / 2, 2.0], [0.4, H / 2, 2.0], [0.0, -H / 2, 2.0]])
    assert np.allclose(structure.thickness(pts), [CAP, SKIN, SKIN])


def test_wall_mid_and_cavity_follow_the_local_thickness(structure):
    # On the upper wall at the cap the mid-surface sits CAP/2 inside the OML.
    at_cap = structure.wall_mid.evaluate(np.array([[0.0, H / 2 - CAP / 2, 1.0]])).phi[0]
    at_skin = structure.wall_mid.evaluate(np.array([[0.4, H / 2 - SKIN / 2, 1.0]])).phi[0]
    assert at_cap == pytest.approx(0.0, abs=1e-12) and at_skin == pytest.approx(0.0, abs=1e-12)
    assert structure.cavity.evaluate(np.array([[0.0, H / 2 - CAP, 1.0]])).phi[0] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("x, y, z, name", [
    (0.4, H / 2 - SKIN / 2, 0.5, "skin"),
    (0.1, H / 2 - CAP / 2, 0.5, "cap_upper"),
    (0.0, 0.0, 0.5, "web"),
    (0.0, H / 2 - CAP / 2, 0.5, "cap_upper"),  # cap outranks web
    (0.35, 0.0, 0.5, "core_aft"),
    (-0.35, -H / 2 + 0.03, 2.0, "stringer"),
    (-0.35, -H / 2 + 0.03, 0.5, None),  # before the stringer starts
    (-0.35, -H / 2 + 0.03, 3.0, "rib_3"),  # the rib outranks everything in its plane
    (-0.2, 0.0, 3.005, "rib_3"),
    (-0.2, 0.0, 3.05, None),
    (0.0, H / 2 + 0.01, 0.5, None),  # outside the OML
])
def test_material_owner(structure, x, y, z, name):
    assert owner(structure, x, y, z) == name


def test_branch_masks(structure):
    web, stringer = structure.branches
    inside = np.array([[0.0, 0.0, 2.0]])
    assert web.mask.evaluate(inside).phi[0] < 0
    near_lower = np.array([[-0.35, -H / 2 + 0.02, 2.0]])
    far = np.array([[-0.35, -H / 2 + 0.08, 2.0]])  # beyond the stringer's depth
    upper = np.array([[-0.35, H / 2 - 0.02, 2.0]])  # wrong side of the clip
    assert stringer.mask.evaluate(near_lower).phi[0] < 0
    assert stringer.mask.evaluate(far).phi[0] > 0
    assert stringer.mask.evaluate(upper).phi[0] > 0
    assert stringer.mask.evaluate(near_lower * [1, 1, 0] + [0, 0, 5.0]).phi[0] > 0


def test_stations(structure):
    assert structure.events == (1.0, 4.0)
    assert structure.transverse_stations == (3.0,)
    assert structure.forced_stations() == (1.0, 3.0, 4.0)
    assert "stringer" in structure.members_at(2.0) and "stringer" not in structure.members_at(5.0)


def test_thickness_can_vary_along_the_span():
    ctx = BuildContext(oml=rect(-0.5, 0.5, -0.1, 0.1), span=(0.0, 2.0))
    skin = WallContribution(ctx.table.add("skin", "skin"), TableLaw([0.0, 2.0], [0.02, 0.01], kind="linear"))
    s = assemble(ctx, [skin])
    pts = np.array([[0.0, 0.1, 0.0], [0.0, 0.1, 1.0], [0.0, 0.1, 2.0]])
    assert np.allclose(s.thickness(pts), [0.02, 0.015, 0.01])


def test_a_new_member_kind_needs_only_a_builder():
    class Doubler(MemberBuilder):
        """Local reinforcement: a partition that reads the cap it sits beside."""

        def build(self, ctx, prior: Sequence[Contribution]):
            assert any(c.spec.name == "skin" for c in prior)
            return [WallContribution(ctx.table.add("doubler", "cap"), ConstantLaw(0.05),
                                     footprint=rect(0.3, 0.4, -1.0, 0.0))]

    ctx = BuildContext(oml=rect(-0.5, 0.5, -0.1, 0.1), span=(0.0, 1.0))
    s = StructureAssembler(BuilderRegistry().register(Skin()).register(Doubler())).assemble(ctx)
    assert s.thickness(np.array([[0.35, -0.1, 0.5]]))[0] == pytest.approx(0.05)


def test_invalid_structures_are_rejected():
    ctx = BuildContext(oml=rect(-0.5, 0.5, -0.1, 0.1), span=(0.0, 1.0))
    a = WallContribution(ctx.table.add("a", "skin"), ConstantLaw(0.01))
    b = WallContribution(ctx.table.add("b", "skin"), ConstantLaw(0.01))
    with pytest.raises(StructureError, match="exactly one base"):
        assemble(ctx, [a, b])
    with pytest.raises(StructureError, match="at least one wall"):
        assemble(ctx, [])
    with pytest.raises(StructureError, match="twice"):
        assemble(ctx, [a, a])
    other = BuildContext(oml=ctx.oml, span=(0.0, 1.0))
    with pytest.raises(StructureError, match="not in the build table"):
        assemble(other, [WallContribution(ctx.table.add("c", "core"), ConstantLaw(0.01))])
    with pytest.raises(ValueError):
        BranchContribution(a.spec, ctx.oml, ctx.oml, ConstantLaw(0.01), free_end=True)
    with pytest.raises(ValueError):
        BuildContext(oml=ctx.oml, span=(1.0, 0.0))


def test_stacked_layer_lies_beneath_the_wall():
    """A cap kept as its own layer under the skin: skin keeps the outer band."""
    ctx = BuildContext(oml=rect(-0.5, 0.5, -0.1, 0.1), span=(0.0, 1.0))
    skin = WallContribution(ctx.table.add("skin", "skin"), ConstantLaw(SKIN))
    cap = WallContribution(ctx.table.add("cap", "cap"), ConstantLaw(CAP),
                           footprint=rect(-0.2, 0.2, 0.0, 1.0), stacked=True)
    s = assemble(ctx, [skin, cap])
    top = 0.1
    assert s.thickness(np.array([[0.0, top, 0.0], [0.4, top, 0.0]])).tolist() == pytest.approx([SKIN + CAP, SKIN])
    assert owner(s, 0.0, top - SKIN / 2, 0.5) == "skin"
    assert owner(s, 0.0, top - SKIN - CAP / 2, 0.5) == "cap"
    assert owner(s, 0.0, top - SKIN - CAP - 0.002, 0.5) is None
    assert owner(s, 0.4, top - SKIN - 0.002, 0.5) is None
    with pytest.raises(ValueError):
        WallContribution(cap.spec, ConstantLaw(CAP), stacked=True)


def test_highest_priority_partition_sets_the_thickness():
    ctx = BuildContext(oml=rect(-0.5, 0.5, -0.1, 0.1), span=(0.0, 1.0))
    skin = WallContribution(ctx.table.add("skin", "skin"), ConstantLaw(0.02))
    thin = WallContribution(ctx.table.add("thin", "cap"), ConstantLaw(0.005), footprint=rect(-0.1, 0.1, 0.0, 1.0))
    s = assemble(ctx, [skin, thin])
    assert s.thickness(np.array([[0.0, 0.1, 0.0]]))[0] == pytest.approx(0.005)


def test_partition_of_equal_priority_replaces_the_base_wall():
    """Two skins of one kind (suction and pressure): each owns its own side."""
    ctx = BuildContext(oml=rect(-0.5, 0.5, -0.1, 0.1), span=(0.0, 1.0))
    upper = WallContribution(ctx.table.add("upper", "skin"), ConstantLaw(0.01))
    lower = WallContribution(ctx.table.add("lower", "skin"), ConstantLaw(0.004),
                             footprint=rect(-1.0, 1.0, -1.0, 0.0))
    s = assemble(ctx, [upper, lower])
    assert owner(s, 0.0, 0.1 - 0.005, 0.5) == "upper"
    assert owner(s, 0.0, -0.1 + 0.002, 0.5) == "lower"
    assert owner(s, 0.0, -0.1 + 0.006, 0.5) is None
