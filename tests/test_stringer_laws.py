"""A stringer's chord fraction and depth can vary along the span."""

import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, Stringer


def _cut(stringer, z):
    foil = Airfoil.naca4("0024")
    stations = [
        Station(0.0, 1.0, foil, pitch_axis=0.3),
        Station(10.0, 1.0, foil, pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-0R", 0.004,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.004,
        web_positions=(0.2, 0.6),
        pivot_fraction=0.3,
        stringers=(stringer,),
    ).section_at(z).stringers[0]


def test_position_and_depth_follow_span_laws():
    stringer = Stringer(
        "aft", 0.4, 0.002, 0.02,
        position_law=lambda z: 0.30 + 0.01 * z,
        depth_law=((0.0, 10.0), (0.02, 0.04)),
    )
    root = _cut(stringer, 0.0)
    tip = _cut(stringer, 10.0)
    assert root.position == pytest.approx(0.30)
    assert tip.position == pytest.approx(0.40)
    assert root.depth == pytest.approx(0.02)
    assert tip.depth == pytest.approx(0.04)


def test_a_fixed_stringer_is_unchanged():
    stringer = Stringer("aft", 0.4, 0.002, 0.02)
    assert _cut(stringer, 4.0).position == pytest.approx(0.4)
    assert _cut(stringer, 4.0).depth == pytest.approx(0.02)
