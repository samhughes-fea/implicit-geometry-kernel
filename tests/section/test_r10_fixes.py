"""Regressions for the feasibility defects found while deriving R10 (SG-038, SG-040)."""

from __future__ import annotations

import pytest

from implicit_geometry_kernel import Airfoil, resolve_system_type
from implicit_geometry_kernel.contract.errors import InfeasibleSectionGeometry
from implicit_geometry_kernel.section.assembly.web_layout import inner_y_at_x, resolve_web_anchors
from implicit_geometry_kernel.section.spec import SectionGeometry, SectionSpec


def _depth(af, t, x):
    try:
        yt, yb = inner_y_at_x(af, t, x, n_samples=20_001, y_search=(-0.05, 0.05))
    except InfeasibleSectionGeometry:
        return 0.0
    return yt - yb


@pytest.mark.parametrize("target", [0.4e-3, 0.6e-3])
def test_shallow_cavity_that_exists_is_kept(target):
    """SG-038: the probe window is padded by the thickness, so a 0.4 mm cavity is found."""
    af = Airfoil.naca4("0006").to_sdf(chord=1.0)
    x = 0.8
    lo, hi = 0.0, 0.5 * _depth(af, 0.0, x)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if _depth(af, mid, x) > target else (lo, mid)
    kept, dropped = resolve_web_anchors(af, lo, [x], return_dropped=True)
    assert len(kept) == 1 and not dropped
    assert kept[0][1] - kept[0][2] == pytest.approx(target, rel=1e-3)


@pytest.mark.parametrize("skin", [0.001, 0.0012, 0.002])
def test_thin_skin_resolves_at_the_floored_skin(skin):
    """SG-040: a skin under 1.25 mm builds; its webs are anchored on the 2.5 mm skin."""
    af = Airfoil.naca4("2412").to_sdf(chord=1.0)
    spec = SectionSpec.from_kwargs(
        resolve_system_type("2D-F-LC-0R"),
        af,
        skin_thickness=skin,
        spar_cap_thickness=0.0,
        shear_web_thickness=0.001,
        pivot_fraction=0.35,
    )
    runner = SectionGeometry(spec).execute()
    walls = runner.get_results().section.wall_recipes
    assert walls["skin_suction_surface"].thickness == 0.0025
    floored = SectionGeometry(
        SectionSpec.from_kwargs(
            resolve_system_type("2D-F-LC-0R"),
            af,
            skin_thickness=0.0025,
            spar_cap_thickness=0.0,
            shear_web_thickness=0.001,
            pivot_fraction=0.35,
        )
    ).execute()
    assert runner.get_layout_resolution().web_y_coords == floored.get_layout_resolution().web_y_coords
