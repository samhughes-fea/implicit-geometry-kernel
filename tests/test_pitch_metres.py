"""Positions in metres from the pitch axis, converted on each station's chord."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, Stringer
from implicit_geometry_kernel.blade.concepts import chord_fraction_from_pitch


AIRFOIL = Airfoil.naca4("0024")


def _stations():
    return [
        Station(0.0, 2.0, AIRFOIL, pitch_axis=0.3),
        Station(10.0, 1.0, AIRFOIL, pitch_axis=0.3),
    ]


def _blade(**kw):
    options = dict(
        pivot_fraction=0.3,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.004,
    )
    options.update(kw)
    return BladeGeometry(_stations(), "2D-F-LC-0R", 0.004, **options)


def test_chord_fraction_uses_this_stations_chord():
    assert chord_fraction_from_pitch(-0.12, 2.0, 0.3) == pytest.approx(0.24)
    assert chord_fraction_from_pitch(-0.12, 1.0, 0.3) == pytest.approx(0.18)


def test_pitch_axis_metres_stay_put_as_the_chord_tapers():
    geo = _blade(web_x_from_pitch=(-0.12, 0.18))
    for z in (0.0, 10.0):
        xs = sorted(float(np.mean(w.xy[:, 0])) for w in geo.section_at(z).walls_of("web"))
        assert xs == pytest.approx([-0.12, 0.18], abs=1e-9)


def test_pitch_axis_metres_follow_a_span_law():
    geo = _blade(web_x_from_pitch=lambda z: (-0.20, 0.05 + 0.01 * z))
    root = sorted(float(np.mean(w.xy[:, 0])) for w in geo.section_at(0.0).walls_of("web"))
    tip = sorted(float(np.mean(w.xy[:, 0])) for w in geo.section_at(10.0).walls_of("web"))
    assert root == pytest.approx([-0.20, 0.05], abs=1e-9)
    assert tip == pytest.approx([-0.20, 0.15], abs=1e-9)


def test_two_foot_webs_from_pitch_metres():
    geo = _blade(web_x_from_pitch=((-0.20, -0.18), (0.10, 0.12)))
    cut = geo.section_at(0.0)
    xs = sorted(float(np.mean(w.xy[:, 0])) for w in cut.walls_of("web"))
    assert xs == pytest.approx([-0.19, 0.11], abs=1e-6)


def test_cap_window_from_pitch_metres():
    geo = _blade(web_positions=(0.2, 0.7), cap_x_from_pitch=((-0.05, 0.10),))
    cap = geo.section_at(10.0).walls_of("cap")[0]
    assert cap.xy[:, 0].min() == pytest.approx(-0.05, abs=1e-6)
    assert cap.xy[:, 0].max() == pytest.approx(0.10, abs=1e-6)


def test_stringer_position_m_uses_the_local_chord():
    stringer = Stringer("aft", 0.5, 0.002, 0.02, position_m=-0.10)
    geo = _blade(web_positions=(0.2, 0.7), stringers=(stringer,))
    assert geo.section_at(0.0).stringers[0].position == pytest.approx(0.3 + -0.10 / 2.0)
    assert geo.section_at(10.0).stringers[0].position == pytest.approx(0.3 + -0.10 / 1.0)


def test_pitch_metres_are_exclusive():
    with pytest.raises(ValueError):
        _blade(web_positions=(0.2, 0.5), web_x_from_pitch=(-0.1, 0.1))
    with pytest.raises(ValueError):
        _blade(cap_windows=((0.2, 0.4),), cap_x_from_pitch=((-0.1, 0.1),))
