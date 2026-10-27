"""Structural families A/B/C/D on MultiCellSection."""

import numpy as np
import pytest

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    SDFGrid,
)
from blade_geometry_model.section_geometry.layout.section_axes import (
    max_thickness_chord_x,
    pitch_axis_x_from_le,
)


@pytest.fixture
def af():
    return AirfoilSDF.from_naca("2412", chord=1.0)


def _area(phi, grid):
    return grid.area(phi)


def test_pitch_and_max_thickness_differ_on_asymmetric(af):
    """2412 is asymmetric enough that pitching (1/3) != max-thickness x."""
    xp = pitch_axis_x_from_le(af, 1.0 / 3.0)
    xm = max_thickness_chord_x(af, n_points=300)
    assert abs(xp - xm) > 1e-3


def test_structural_D_matches_twin_web_default(af):
    m_default = MultiCellSection.twin_web(af, twist_angle=0.0)
    m_d = MultiCellSection(
        af,
        web_x_positions=[0.2, 0.5],
        structural_family="D",
        twist_angle=0.0,
    )
    assert m_default.structural_family == "D"
    assert m_d.structural_family == "D"
    grid = SDFGrid.from_airfoil(af, nx=180, ny=96)
    a0 = _area(grid.eval(m_default["spar_cap_suction_surface"]), grid)
    a1 = _area(grid.eval(m_d["spar_cap_suction_surface"]), grid)
    assert abs(a0 - a1) / max(a0, 1e-12) < 0.02


def test_structural_A_has_no_spar_caps(af):
    m = MultiCellSection(
        af,
        web_x_positions=[0.33, 0.66],
        structural_family="A",
        twist_angle=0.0,
    )
    assert "spar_cap_suction_surface" not in m.labels
    assert "spar_cap_pressure_surface" not in m.labels
    assert "shear_web_0" in m.labels


def test_structural_B_pitching_vs_max_thickness(af):
    m_pitch = MultiCellSection(
        af,
        web_x_positions=[0.25, 0.55],
        structural_family="B",
        fixed_cap_anchor="pitching",
        pitch_fraction_of_chord_from_le=1.0 / 3.0,
        fixed_cap_chord_half_width=0.06,
        twist_angle=0.0,
    )
    m_max = MultiCellSection(
        af,
        web_x_positions=[0.25, 0.55],
        structural_family="B",
        fixed_cap_anchor="max_thickness",
        fixed_cap_chord_half_width=0.06,
        twist_angle=0.0,
    )
    grid = SDFGrid.from_airfoil(af, nx=200, ny=100)
    # Centres differ → upper cap areas should not match exactly
    ap = _area(grid.eval(m_pitch["spar_cap_suction_surface"]), grid)
    am = _area(grid.eval(m_max["spar_cap_suction_surface"]), grid)
    assert ap > 0 and am > 0
    assert abs(ap - am) / max(ap, 1e-12) > 0.01


def test_structural_C_discrete_caps(af):
    m = MultiCellSection(
        af,
        web_x_positions=[0.3, 0.55, 0.78],
        structural_family="C",
        discrete_cap_chord_half_width=0.035,
        twist_angle=0.0,
    )
    grid = SDFGrid.from_airfoil(af, nx=220, ny=110)
    a = _area(grid.eval(m["spar_cap_suction_surface"]), grid)
    assert a > 0


@pytest.mark.parametrize("web_xs", [
    [0.45],
    [0.25, 0.45, 0.65],
])
def test_structural_C_cap_recipes_centered_on_webs(af, web_xs):
    """SG-013: each CapWall midpoint stays on its web station."""
    m = MultiCellSection(
        af,
        web_x_positions=web_xs,
        structural_family="C",
        discrete_cap_chord_half_width=0.05,
        skin_thickness=0.006,
        web_thickness=0.008,
        cap_height=0.014,
        twist_angle=0.0,
    )
    recipes = m.spar_cap_wall_recipes
    for i, web_x in enumerate(web_xs):
        key = f"spar_cap_suction_surface_{i}"
        assert key in recipes, sorted(recipes)
        wall = recipes[key]
        mid = 0.5 * (float(wall.x_start) + float(wall.x_end))
        assert mid == pytest.approx(float(web_x), abs=2e-3)
        assert (web_x - wall.x_start) == pytest.approx(
            wall.x_end - web_x, abs=2e-3
        )


def test_B_C_require_airfoil_sdf_not_generic_callable():
    class Dummy:
        chord = 1.0

        def __call__(self, x, y):
            return np.ones(np.shape(x)) * 1.0

    with pytest.raises(TypeError, match="AirfoilSDF"):
        MultiCellSection(
            Dummy(),
            web_x_positions=[0.4],
            structural_family="B",
        )


def test_zero_web_family_A_builds_skin_and_one_core(af):
    m = MultiCellSection(af, web_x_positions=[], structural_family="A")
    assert m.n_webs == 0
    assert m.n_cells == 1
    assert "spar_cap_suction_surface" not in m.labels
    assert "spar_cap_pressure_surface" not in m.labels
    assert "sandwich_core_0" in m.labels
    grid = SDFGrid.from_airfoil(af, nx=180, ny=96)
    assert _area(grid.eval(m["sandwich_core_0"]), grid) > 0


def test_zero_web_family_B_builds_caps_and_one_core(af):
    m = MultiCellSection(
        af,
        web_x_positions=[],
        structural_family="B",
        fixed_cap_chord_half_width=0.06,
    )
    assert m.n_webs == 0
    assert m.n_cells == 1
    assert "spar_cap_suction_surface" in m.labels
    assert "spar_cap_pressure_surface" in m.labels
    assert "sandwich_core_0" in m.labels
    grid = SDFGrid.from_airfoil(af, nx=200, ny=100)
    assert _area(grid.eval(m["spar_cap_suction_surface"]), grid) > 0
    assert _area(grid.eval(m["sandwich_core_0"]), grid) > 0


@pytest.mark.parametrize("family", ["C", "D"])
def test_zero_web_family_C_D_rejected(af, family):
    with pytest.raises(ValueError, match="empty"):
        MultiCellSection(af, web_x_positions=[], structural_family=family)
