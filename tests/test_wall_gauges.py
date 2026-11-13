"""Cap sides and webs can carry their own thickness."""

import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station


def _cut(**kw):
    foil = Airfoil.naca4("0024")
    stations = [
        Station(0.0, 1.0, foil, pitch_axis=0.3),
        Station(4.0, 1.0, foil, pitch_axis=0.3),
    ]
    options = dict(
        pivot_fraction=0.3,
        spar_cap_thickness=(0.008, 0.014),
        shear_web_thickness=(0.005, 0.009),
        web_positions=(0.2, 0.55),
    )
    options.update(kw)
    return BladeGeometry(stations, "2D-F-LC-0R", 0.004, **options).section_at(1.0)


def test_each_cap_side_and_each_web_keeps_its_thickness():
    cut = _cut()
    caps = {w.label: w.thickness for w in cut.walls_of("cap")}
    assert caps["spar_cap_suction_surface"] == pytest.approx(0.008)
    assert caps["spar_cap_pressure_surface"] == pytest.approx(0.014)
    webs = {w.label: w.thickness for w in cut.walls_of("web")}
    assert webs["shear_web_0"] == pytest.approx(0.005)
    assert webs["shear_web_1"] == pytest.approx(0.009)


def test_a_span_law_can_return_a_gauge_pair():
    cut = _cut(
        spar_cap_thickness=lambda z: (0.007, 0.011),
        shear_web_thickness=lambda z: (0.006, 0.008),
    )
    caps = {w.label: w.thickness for w in cut.walls_of("cap")}
    assert caps["spar_cap_suction_surface"] == pytest.approx(0.007)
    assert caps["spar_cap_pressure_surface"] == pytest.approx(0.011)
