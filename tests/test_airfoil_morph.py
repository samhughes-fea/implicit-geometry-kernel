"""The station-table airfoil morph is an explicit option (G7)."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station


def _blade(rule):
    ellipse, naca = Airfoil.ellipse(1.0), Airfoil.naca4("0012")
    stations = [
        Station(0.0, 2.0, ellipse, pitch_axis=0.5),
        Station(10.0, 1.0, naca, pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "0A-LC-0R", 0.004, pivot_fraction=0.3, airfoil_morph=rule,
    )


def test_auto_keeps_the_smoothstep_on_a_circle_root():
    geo = _blade("auto")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, geo.n_airfoil_points)))
    *_, (_x, upper, _lower) = geo.interpolate(2.5)
    w = 0.25
    s = w * w * (3.0 - 2.0 * w)
    ellipse, naca = geo.stations[0].airfoil, geo.stations[1].airfoil
    smooth = (1.0 - s) * ellipse.surface("upper", x) + s * naca.surface("upper", x)
    linear = (1.0 - w) * ellipse.surface("upper", x) + w * naca.surface("upper", x)
    assert np.max(np.abs(upper - smooth)) < np.max(np.abs(upper - linear))


def test_linear_option_does_not_smoothstep():
    geo = _blade("linear")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, geo.n_airfoil_points)))
    *_, (_x, upper, _lower) = geo.interpolate(2.5)
    w = 0.25
    ellipse, naca = geo.stations[0].airfoil, geo.stations[1].airfoil
    linear = (1.0 - w) * ellipse.surface("upper", x) + w * naca.surface("upper", x)
    # the trailing-edge gap is still the linear gap; the rest of the surface matches
    assert np.max(np.abs(upper[:-1] - linear[:-1])) < 1e-12


def test_smoothstep_option_on_two_foils():
    a, b = Airfoil.naca4("0018"), Airfoil.naca4("0012")
    stations = [Station(0.0, 1.0, a, pitch_axis=0.3), Station(4.0, 1.0, b, pitch_axis=0.3)]
    geo = BladeGeometry(stations, "0A-LC-0R", 0.003, pivot_fraction=0.3, airfoil_morph="smoothstep")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, geo.n_airfoil_points)))
    *_, (_x, upper, _lower) = geo.interpolate(1.0)
    w = 0.25
    s = w * w * (3.0 - 2.0 * w)
    expect = (1.0 - s) * a.surface("upper", x) + s * b.surface("upper", x)
    assert np.max(np.abs(upper[:-1] - expect[:-1])) < 1e-12


def test_unknown_morph_is_rejected():
    with pytest.raises(ValueError, match="airfoil_morph"):
        _blade("cubic")
