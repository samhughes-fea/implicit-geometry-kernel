"""A chord-fraction web can tilt away from the thickness direction."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station


def _blade(**kw):
    foil = Airfoil.naca4("0024")
    stations = [
        Station(0.0, 1.0, foil, pitch_axis=0.3),
        Station(4.0, 1.0, foil, pitch_axis=0.3),
    ]
    options = dict(
        pivot_fraction=0.3,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.004,
        web_positions=(0.25, 0.55),
    )
    options.update(kw)
    return BladeGeometry(stations, "2D-CN-LC-0R", 0.004, **options)


def _angle(web):
    """Degrees from the thickness direction toward the trailing edge, suction end up."""
    d = np.asarray(web.xy[-1] - web.xy[0], float)
    if web.xy[0, 1] < web.xy[-1, 1]:
        d = -d
    return float(np.degrees(np.arctan2(d[0], abs(d[1]))))


def test_zero_tilt_stays_chord_normal():
    plain = _blade()
    tilted = _blade(web_tilt_deg=0.0)
    for a, b in zip(plain.section_at(1.0).walls_of("web"), tilted.section_at(1.0).walls_of("web")):
        assert _angle(a) == pytest.approx(_angle(b), abs=1e-6)


def test_tilt_rotates_the_web_toward_the_trailing_edge():
    webs = _blade(web_tilt_deg=12.0).section_at(1.0).walls_of("web")
    for web in webs:
        assert _angle(web) == pytest.approx(12.0, abs=0.5)


def test_each_web_has_its_own_tilt():
    webs = _blade(web_tilt_deg=(0.0, 15.0)).section_at(1.0).walls_of("web")
    angles = sorted(abs(_angle(w)) for w in webs)
    assert angles[0] == pytest.approx(0.0, abs=0.5)
    assert angles[1] == pytest.approx(15.0, abs=0.5)


def test_two_foot_webs_reject_a_tilt():
    with pytest.raises(ValueError):
        _blade(
            web_positions=lambda z: ((0.2, 0.22), (0.5, 0.55)),
            web_tilt_deg=10.0,
        ).section_at(1.0)
