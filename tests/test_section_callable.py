"""Callable section(z) lofted to a tolerance, without a second airfoil morph."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, SectionSample, Station
from implicit_geometry_kernel.blade.loft.section_callable import placed_oml, polyline_gap


def _sample(z, chord, twist, pivot, offset, foil, x):
    return SectionSample(
        float(z), float(chord), float(twist), float(pivot),
        (float(offset[0]), float(offset[1])),
        np.asarray(x, float),
        np.asarray(foil.surface("upper", x), float),
        np.asarray(foil.surface("lower", x), float),
    )


def test_linear_chord_needs_no_extra_stations():
    foil = Airfoil.naca4("0012")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 81)))

    def section(z):
        return _sample(z, 1.0 + 0.1 * z, 2.0 * z, 0.3, (0.01 * z, -0.02 * z), foil, x)

    geo = BladeGeometry.from_section(section, (0.0, 10.0), 1e-4, "0A-LC-0R", 0.002, n_airfoil_points=81)
    assert [s.z for s in geo.stations] == pytest.approx([0.0, 10.0])
    assert geo._section_samples is not None


def test_quadratic_chord_is_refined_to_the_tolerance():
    foil = Airfoil.naca4("0018")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 61)))

    def section(z):
        return _sample(z, 1.0 + z * z, 8.0 * np.sin(0.5 * z), 0.25 + 0.02 * z, (0.0, 0.05 * z), foil, x)

    tol = 1e-4
    geo = BladeGeometry.from_section(section, (0.0, 2.0), tol, "0A-LC-0R", 0.002, n_airfoil_points=61)
    assert len(geo.stations) > 2
    for z in np.linspace(0.0, 2.0, 17):
        gap = polyline_gap(placed_oml(section(z)), geo.oml_polyline(float(z))[:, :2])
        assert gap <= tol + 1e-12


def test_callable_path_does_not_morph_again():
    circle, foil = Airfoil.ellipse(1.0), Airfoil.naca4("0012")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 49)))
    circle_u = np.asarray(circle.surface("upper", x), float)
    circle_l = np.asarray(circle.surface("lower", x), float)
    foil_u = np.asarray(foil.surface("upper", x), float)
    foil_l = np.asarray(foil.surface("lower", x), float)

    def section(z):
        # The consumer's blend, already evaluated. Not Airfoil.blend, which refits a spline.
        w = float(np.clip(z / 4.0, 0.0, 1.0))
        s = w * w * (3.0 - 2.0 * w)
        return SectionSample(
            float(z), 2.0 - 0.2 * z, 0.0, 0.5, (0.0, 0.0), x,
            (1.0 - s) * circle_u + s * foil_u,
            (1.0 - s) * circle_l + s * foil_l,
        )

    def fail_morph(*_args, **_kwargs):
        raise AssertionError("_morph_surfaces must not run on a callable section")

    import implicit_geometry_kernel.blade.loft.loft_cut as loft_cut

    original = loft_cut._morph_surfaces
    loft_cut._morph_surfaces = fail_morph
    try:
        geo = BladeGeometry.from_section(
            section, (0.0, 4.0), 1e-4, "0A-LC-0R", 0.002, knots=(1.0,), n_airfoil_points=49,
        )
        assert any(abs(s.z - 1.0) < 1e-9 for s in geo.stations)
        for z in (0.0, 0.3, 1.0, 2.2, 4.0):
            gap = polyline_gap(placed_oml(section(z)), geo.oml_polyline(z)[:, :2])
            assert gap <= 1e-4 + 1e-12
        # the station table still morphs
        loft_cut._morph_surfaces = original
        stations = [
            Station(0.0, 2.0, circle, pitch_axis=0.5),
            Station(4.0, 1.0, foil, pitch_axis=0.5),
        ]
        table = BladeGeometry(stations, "0A-LC-0R", 0.002, pivot_fraction=0.5, n_airfoil_points=49)
        assert table._section_samples is None
        called = {"n": 0}

        def count(*args, **kwargs):
            called["n"] += 1
            return original(*args, **kwargs)

        loft_cut._morph_surfaces = count
        table.interpolate(1.0)
        assert called["n"] == 1
    finally:
        loft_cut._morph_surfaces = original


def test_tuple_return_and_endpoint_exact():
    foil = Airfoil.naca4("2412")

    def section(z):
        return (1.5, 3.0, 0.3, (0.1, -0.2), foil)

    geo = BladeGeometry.from_section(section, (0.0, 1.0), 1e-3, "0A-LC-0R", 0.002)
    chord, twist, pitch, offset, _surfaces = geo.interpolate(0.0)
    assert (chord, twist, pitch) == pytest.approx((1.5, 3.0, 0.3))
    assert offset == pytest.approx((0.1, -0.2))
    chord, twist, pitch, offset, _surfaces = geo.interpolate(1.0)
    assert (chord, twist, pitch) == pytest.approx((1.5, 3.0, 0.3))
