"""Names and geometry that downstream code written against 0.2.x relies on."""

import pytest

import implicit_geometry_kernel as k
from implicit_geometry_kernel import Airfoil


def test_shear_web_names():
    from implicit_geometry_kernel.section.walls.shear_web import ShearWebComponent

    assert k.ShearWeb is not ShearWebComponent


@pytest.mark.parametrize("code", ["0015", "0021"])
def test_naca4_n_is_points_per_surface(code):
    # ``n`` is points per surface, as in 0.2.0. At half that density the spline through
    # the points overshoots the nose and the leading edge opens by about 0.05% chord.
    a = Airfoil.naca4(code)
    gap = a.surface("upper", 0.0) - a.surface("lower", 0.0)
    assert abs(gap) < 2.5e-4


# 0.2.0 geometry of a NACA 0024 two-web blade cut at z = 3 m, recorded from v0.2.0
# with BladeGeometry(stations, "2D-F", 0.003, 0.008, 0.002, web_positions=(0.18, 0.48)).
_V020_WALL_LENGTHS = {
    "shear_web_0": 0.221203989,
    "shear_web_1": 0.213341770,
    "spar_cap_suction_surface": 0.300867922,
}
_V020_CELL_AREAS = (0.029221933, 0.069006524, 0.062002916)


def _migrated_cut():
    from implicit_geometry_kernel import BladeGeometry, Station

    af = Airfoil.naca4("0024")
    st = [Station(0.0, 1.0, af, pitch_axis=0.3), Station(6.0, 1.0, af, pitch_axis=0.3)]
    geo = BladeGeometry(st, "2D-F-LC-0R", 0.003, 0.008, 0.002, web_positions=(0.18, 0.48), pivot_fraction=0.3)
    return geo.section_at(3.0)


def _length(points):
    import numpy as np

    p = np.asarray(points, float)[:, :2]
    return float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())


def _area(points):
    import numpy as np

    p = np.asarray(points, float)[:, :2]
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))))


def test_migrated_call_reproduces_0_2_geometry():
    # docs/MIGRATING.md: full layout token plus pivot_fraction gives 0.2.0 geometry
    cut = _migrated_cut()
    walls = {w.label: w for w in cut.walls}
    for label, length in _V020_WALL_LENGTHS.items():
        assert _length(walls[label].points) == pytest.approx(length, rel=5e-4)
    areas = sorted(_area(c.points) for c in cut.cells)
    assert areas == pytest.approx(sorted(_V020_CELL_AREAS), rel=5e-4)


def test_short_layout_code_names_the_full_token():
    from implicit_geometry_kernel import resolve_system_type

    with pytest.raises(KeyError, match="2D-F-LC-0R"):
        resolve_system_type("2D-F")


def test_blade_geometry_requires_pivot_fraction():
    from implicit_geometry_kernel import BladeGeometry, Station

    af = Airfoil.naca4("0024")
    st = [Station(0.0, 1.0, af), Station(6.0, 1.0, af)]
    with pytest.raises(TypeError, match="pivot_fraction"):
        BladeGeometry(st, "2D-F-LC-0R", 0.003, 0.008, 0.002)


def test_system_type_keys_at_top_level():
    import implicit_geometry_kernel as kernel

    assert "2D-F-LC-0R" in kernel.SYSTEM_TYPE_KEYS
