"""SG-033: a flapwise web at 90 degrees of twist runs along the chord."""

from __future__ import annotations

import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.contract import InfeasibleSectionGeometry


def _blade(layout: str, twist: float) -> BladeGeometry:
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0018"), twist=twist, pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), twist=twist, pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations,
        layout,
        skin_thickness=0.004,
        pivot_fraction=0.3,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.004,
        web_positions=(0.18, 0.48),
    )


@pytest.mark.parametrize("twist", [90.0, -90.0, 270.0, 89.5])
def test_flapwise_webs_along_the_chord_raise_a_clear_error(twist: float) -> None:
    with pytest.raises(InfeasibleSectionGeometry, match="chord-normal") as caught:
        _blade("2D-F-LC-0R", twist).section_at(4.0)
    assert caught.value.web_index == 0


@pytest.mark.parametrize("twist", [60.0, 88.9])
def test_flapwise_webs_short_of_the_limit_still_build(twist: float) -> None:
    section = _blade("2D-F-LC-0R", twist).section_at(4.0)
    assert section.station_webs == 2


def test_chord_normal_webs_build_at_ninety_degrees() -> None:
    section = _blade("2D-CN-LC-0R", 90.0).section_at(4.0)
    assert section.station_webs == 2
