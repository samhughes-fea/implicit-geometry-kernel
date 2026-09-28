"""WP3: anchored wall geometry, collision checks, and strip placement stamps."""

from __future__ import annotations

import numpy as np
import pytest

from blade_structure_utilities.geometry_utilities.section_contract.shell_midline_strip import (
    ShellMidlineStrip,
)
from blade_structure_utilities.geometry_utilities.section_geometry.placement import (
    PlacementGeometryError,
    anchored_wall_geometry,
    anchored_wall_geometry_from_policy,
    assert_anchored_wall_feasible,
    find_placement_collisions,
    midsurface_offset_from_reference_m,
    require_placement_stamp,
    stamp_strip_placement_fields,
)
from blade_structure_utilities.geometry_utilities.section_geometry.placement.anchoring import (
    airfoil_inward_depth_for_midsurface,
)
from blade_structure_utilities.geometry_utilities.section_geometry.recipes.sdf_ir.outer_skin import (
    OuterSkin,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.airfoil import AirfoilSDF
from blade_structure_utilities.material_utilities.layup import (
    DEFAULT_PLY_PLACEMENT,
    SKIN_CENTRED,
    SKIN_FROM_INNER_MOULD_LINE,
    SKIN_FROM_OUTER_MOULD_LINE,
    SPAR_CAP_FROM_SKIN_BOND_LINE,
    WEB_FROM_MIDSURFACE,
)


@pytest.mark.parametrize(
    "policy, mid, free, second",
    [
        (SKIN_FROM_OUTER_MOULD_LINE, 0.002, 0.004, None),
        (SKIN_FROM_INNER_MOULD_LINE, 0.002, 0.004, None),
        (SPAR_CAP_FROM_SKIN_BOND_LINE, 0.006, 0.012, None),
        (WEB_FROM_MIDSURFACE, 0.0, 0.003, -0.003),
        (SKIN_CENTRED, 0.0, 0.002, -0.002),
    ],
)
def test_anchored_offsets_match_policy(policy, mid, free, second):
    thickness = abs(free) * 2 if second is not None else free
    wall = anchored_wall_geometry_from_policy(policy, thickness)
    assert wall.midsurface_offset_m == pytest.approx(mid)
    assert wall.free_face_offset_m == pytest.approx(free)
    assert wall.second_free_face_offset_m == (
        None if second is None else pytest.approx(second)
    )
    assert_anchored_wall_feasible(wall)


def test_default_policies_match_deck_table():
    assert DEFAULT_PLY_PLACEMENT["skin"] is SKIN_FROM_OUTER_MOULD_LINE
    assert DEFAULT_PLY_PLACEMENT["spar"] is SPAR_CAP_FROM_SKIN_BOND_LINE
    assert DEFAULT_PLY_PLACEMENT["web"] is WEB_FROM_MIDSURFACE


def test_collision_raises_on_negative_cavity():
    skin = anchored_wall_geometry("outer_mould_line", "inward", 0.004)
    hits = find_placement_collisions(
        [skin], skin_cavity_clearance_m=0.02, cavity_depth_m=0.01
    )
    assert any(h.code == "negative_cavity" for h in hits)


def test_collision_raises_on_cap_intrusion():
    cap = anchored_wall_geometry("skin_bond_line", "into_cavity", 0.02)
    hits = find_placement_collisions(
        [cap], cap_into_cavity_m=0.03, cavity_depth_m=0.01
    )
    assert any(h.code == "cap_skin_separation" for h in hits)


def test_require_placement_stamp_refuses_partial_claim():
    with pytest.raises(ValueError, match="placement_digest"):
        require_placement_stamp(
            placement_policy_id="skin_from_outer_mould_line",
            placement_digest=None,
            reference_face="outer_mould_line",
        )


def test_stamp_strip_placement_writes_all_fields():
    strip = ShellMidlineStrip(
        label="skin_suction_surface",
        kind="skin",
        midline_twisted_pitch_axis_blade_frame=np.array([[0.0, 0.0], [1.0, 0.0]]),
        thickness_m=0.004,
    )
    stamped = stamp_strip_placement_fields(
        strip,
        placement_policy_id=SKIN_FROM_OUTER_MOULD_LINE.name,
        placement_digest=SKIN_FROM_OUTER_MOULD_LINE.digest,
        reference_face=SKIN_FROM_OUTER_MOULD_LINE.reference_face.value,
    )
    assert stamped.placement_policy_id == SKIN_FROM_OUTER_MOULD_LINE.name
    assert stamped.placement_digest == SKIN_FROM_OUTER_MOULD_LINE.digest
    assert stamped.reference_face == "outer_mould_line"


def test_fixed_oml_midline_moves_inward_with_thickness():
    """Fixed-OML: authoring contour (OML) unchanged; mid depth grows with t."""
    af = AirfoilSDF.from_naca("0012", n_points=120, chord=1.0)
    thin = OuterSkin(af, 0.002)
    thick = OuterSkin(af, 0.006)
    assert thin.reference_face == "outer_mould_line"
    assert thick._mid_depth == pytest.approx(0.003)
    assert thin._mid_depth == pytest.approx(0.001)
    xs = np.linspace(0.1, 0.8, 5)
    for x in xs:
        assert thin.outer_sdf(x, 0.0) == pytest.approx(thick.outer_sdf(x, 0.0), abs=1e-9)


def test_fixed_iml_mid_depth_moves_toward_fluid_when_thickness_grows():
    """Fixed IML at authoring depth 0.006: thicker wall → shallower mid from OML."""
    t0 = 0.006
    thin = airfoil_inward_depth_for_midsurface(
        "inner_mould_line", 0.004, reference_iml_depth_m=t0
    )
    thick = airfoil_inward_depth_for_midsurface(
        "inner_mould_line", 0.006, reference_iml_depth_m=t0
    )
    assert thin == pytest.approx(t0 - 0.002)
    assert thick == pytest.approx(t0 - 0.003)
    assert thin > thick


def test_web_midsurface_offset_is_zero():
    assert midsurface_offset_from_reference_m("midsurface", "outward_symmetric", 0.008) == 0.0


def test_invalid_thickness_raises():
    with pytest.raises(ValueError):
        anchored_wall_geometry("outer_mould_line", "inward", 0.0)
    from blade_structure_utilities.geometry_utilities.section_geometry.placement.anchoring import (
        AnchoredWallGeometry,
    )

    bad = AnchoredWallGeometry(
        reference_face="outer_mould_line",
        build_direction="inward",
        thickness_m=0.001,
        reference_offset_m=0.0,
        free_face_offset_m=-0.001,
        midsurface_offset_m=0.0005,
    )
    with pytest.raises(PlacementGeometryError):
        assert_anchored_wall_feasible(bad)
