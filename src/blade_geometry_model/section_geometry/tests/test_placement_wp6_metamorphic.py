"""WP6 metamorphic / geometry integration for surface-aware placement."""

from __future__ import annotations

import numpy as np
import pytest

from blade_structure_utilities.geometry_utilities.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.placement import (
    anchored_wall_geometry,
    default_construction_for_kind,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.airfoil import AirfoilSDF
from blade_structure_utilities.material_utilities.layup import (
    SKIN_FROM_OUTER_MOULD_LINE,
    stack_for_placement,
)


@pytest.mark.parametrize("layout_key", ["2D-F", "3D-F", "4D-F"])
def test_default_skin_cap_web_midlines_exist(layout_key: str) -> None:
    af = AirfoilSDF.from_naca("0012", n_points=160, chord=1.0)
    section = build_section_view(
        af,
        resolve_system_type(layout_key),
        skin_thickness=0.003,
        web_thickness=0.003,
        cap_height=0.012,
    )
    comps = dict(getattr(section, "components_unrotated", {}) or {})
    caps = dict(getattr(section, "spar_cap_components_unrotated", {}) or {})
    assert any("skin" in str(k).lower() or "outer" in str(k).lower() for k in comps)
    assert caps, f"{layout_key} produced no caps"
    assert any("web" in str(k).lower() for k in comps)


def test_reversing_stack_preserves_mirror_symmetry():
    """Metamorphic: reversing a placement stack is the same laminate (B=0)."""
    unit = (45.0, -45.0)
    stack = stack_for_placement(unit, 3, SKIN_FROM_OUTER_MOULD_LINE)
    assert stack == tuple(reversed(stack))


def test_reference_relabel_preserves_offset_magnitudes():
    """Equivalent constructions keep |free-face| = thickness for one-sided walls."""
    a = anchored_wall_geometry("outer_mould_line", "inward", 0.005)
    b = anchored_wall_geometry("skin_bond_line", "into_cavity", 0.005)
    assert abs(a.free_face_offset_m) == pytest.approx(abs(b.free_face_offset_m))
    assert a.midsurface_offset_m == pytest.approx(b.midsurface_offset_m)


def test_default_construction_tokens_match_strip_kinds():
    for kind in ("skin", "cap", "web"):
        ref, direction, policy_id = default_construction_for_kind(kind)
        assert ref and direction and policy_id
