"""Station twist is nose up by default; nose_down is the other sense."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station


def _stations(twist):
    foil = Airfoil.naca4("0012")
    return [
        Station(0.0, 1.0, foil, twist=twist, pitch_axis=0.3),
        Station(4.0, 1.0, foil, twist=twist, pitch_axis=0.3),
    ]


def _cut(twist, sense="nose_up"):
    return BladeGeometry(
        _stations(twist), "2D-F-LC-0R", 0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
        pivot_fraction=0.3, twist_sense=sense,
    ).section_at(1.0)


def test_nose_up_is_the_default_and_lifts_the_leading_edge():
    geo = BladeGeometry(
        _stations(15.0), "2D-F-LC-0R", 0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3,
    )
    assert geo.twist_sense == "nose_up"
    le = _cut(15.0).to_blade([-0.3], [0.0])[0]
    assert le[1] > 0.0


def test_nose_down_places_the_opposite_angle():
    down = _cut(15.0, "nose_down")
    opposite = _cut(-15.0)
    assert down.twist == pytest.approx(-15.0)
    assert down.to_blade([-0.3], [0.0])[0, 1] < 0.0
    for wall in down.walls:
        other = opposite.wall(wall.label)
        assert np.allclose(wall.xy, other.xy, atol=1e-9)


def test_unknown_twist_sense_is_rejected():
    with pytest.raises(ValueError):
        BladeGeometry(
            _stations(0.0), "2D-F-LC-0R", 0.004,
            spar_cap_thickness=0.008, shear_web_thickness=0.004,
            pivot_fraction=0.3, twist_sense="left",
        )
