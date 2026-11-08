"""Named shell-vertex locators along span."""

import numpy as np

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.extract.curves import feature_curve, locate_shell_vertex


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)


def test_locate_shell_vertex_and_feature_curve_agree():
    geo = _geo()
    zs = np.linspace(0.0, 10.0, 7)
    vid = ("TE", "suction")
    a = locate_shell_vertex(geo, vid, zs)
    b = feature_curve(geo, vid, zs)
    assert len(a) == len(zs)
    np.testing.assert_allclose(a, b, atol=1e-9)
    assert np.all(np.diff(a[:, 2]) > 0.0)


def test_web_foot_stays_on_section():
    geo = _geo()
    zs = np.linspace(0.0, 10.0, 5)
    pts = locate_shell_vertex(geo, ("anchor", "web_0", "suction"), zs)
    assert len(pts) == 5
    for p, z in zip(pts, zs):
        cut = geo.section_at(float(z))
        xy = np.vstack([w.xy for w in cut.walls_of("web")])
        dist = float(np.min(np.linalg.norm(xy - p[:2], axis=1)))
        assert dist < 0.05 * cut.chord
