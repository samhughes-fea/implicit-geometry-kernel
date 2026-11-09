"""PA regions and hard-coded LC/LE/PA positions."""

from __future__ import annotations

import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.concepts.web_anchoring import (
    lc_web_fractions,
    le_web_offsets_m,
    pa_regions,
    pa_web_x_from_pitch,
)


def _stations():
    foil = Airfoil.naca4("0012")
    return [
        Station(0.0, 1.0, foil, twist=0.0, pitch_axis=0.5),
        Station(5.0, 2.0, foil, twist=0.0, pitch_axis=0.5),
        Station(10.0, 0.6, foil, twist=0.0, pitch_axis=0.5),
    ]


def test_lc_and_le_are_authored_lists():
    fracs = (0.15, 0.50)
    assert lc_web_fractions(fracs) == pytest.approx(fracs)
    assert le_web_offsets_m(fracs, 2.0) == pytest.approx((0.30, 1.00))


def test_pa_regions_are_root_band_and_te_shelf():
    st = _stations()
    pitch = lambda z: 0.5
    r = pa_regions(st, pitch)
    assert r.primary == pytest.approx((-0.5, 0.5))
    assert r.secondary == pytest.approx((0.5, 1.0))
    assert r.z_max_chord == pytest.approx(5.0)


def test_two_web_pa_straddles_pitch_axis_inside_root():
    st = _stations()
    xs = pa_web_x_from_pitch(st, lambda z: 0.5, 2)
    assert xs == pytest.approx((-0.5 + 1.0 / 3.0, -0.5 + 2.0 / 3.0))
    r = pa_regions(st, lambda z: 0.5)
    assert r.primary[0] < xs[0] < 0.0
    assert 0.0 < xs[1] < r.primary[1]


def test_three_web_pa_is_uniform_root_interiors():
    st = _stations()
    xs = pa_web_x_from_pitch(st, lambda z: 0.5, 3)
    assert xs == pytest.approx((-0.25, 0.0, 0.25))


def test_four_web_pa_last_slot_is_secondary():
    st = _stations()
    xs = pa_web_x_from_pitch(st, lambda z: 0.5, 4)
    r = pa_regions(st, lambda z: 0.5)
    assert xs == pytest.approx((-0.25, 0.0, 0.25, 0.5))
    assert all(r.primary[0] <= x <= r.primary[1] for x in xs[:3])
    dx = (r.primary[1] - r.primary[0]) / 4
    assert xs[1] - xs[0] == pytest.approx(dx)
    assert xs[2] - xs[1] == pytest.approx(dx)
    assert xs[3] - xs[2] == pytest.approx(dx)
    assert xs[3] == pytest.approx(r.primary[1])


def test_five_web_pa_second_secondary_is_two_dx():
    st = _stations()
    xs = pa_web_x_from_pitch(st, lambda z: 0.5, 5)
    r = pa_regions(st, lambda z: 0.5)
    dx = (r.primary[1] - r.primary[0]) / 4
    assert xs == pytest.approx((-0.25, 0.0, 0.25, 0.5, 0.75))
    assert xs[4] - xs[2] == pytest.approx(2 * dx)
    assert xs[4] == pytest.approx(r.primary[1] + dx)


def test_blade_geometry_uses_pa_helper():
    st = _stations()
    geo = BladeGeometry(
        st, "2D-F-PA-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
        pivot_fraction=0.5,
    )
    assert geo._web_x_from_pitch == pytest.approx(pa_web_x_from_pitch(st, lambda z: 0.5, 2))
