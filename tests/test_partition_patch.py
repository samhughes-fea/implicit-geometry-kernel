"""Smooth partition patches with analytic derivatives and curvature (decision 0017)."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.patch import PartitionPatch
from implicit_geometry_kernel.structure import UnderSampledSpanLaw

SKIN_LE = "skin_suction_surface:leading_edge"


def _cylinder():
    foil = Airfoil.ellipse(1.0)
    stations = [Station(0.0, 1.0, foil, twist=0.0, pitch_axis=0.5), Station(4.0, 1.0, foil, twist=0.0, pitch_axis=0.5)]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                         shear_web_thickness=0.004, pivot_fraction=0.5)


def _twisted(n):
    foil = Airfoil.naca4("0018")
    zs = np.linspace(0.0, 6.0, n)
    stations = [Station(float(z), float(1.0 - 0.06 * z), foil, twist=float(12.0 * np.sin(z / 2.0)), pitch_axis=0.3)
                for z in zs]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                         shear_web_thickness=0.004, pivot_fraction=0.3, span_interp="pchip")


@pytest.fixture(scope="module")
def cylinder_skin():
    return PartitionPatch(_cylinder(), SKIN_LE, np.linspace(0.0, 4.0, 9))


def test_a_cylinder_skin_is_singly_curved_with_the_mid_surface_radius(cylinder_skin):
    radius = 0.5 - 0.002
    for t in (0.1, 0.3, 0.5, 0.7, 0.8):
        sp = cylinder_skin.at(t, 2.0)
        assert sp.kind == "SINGLY"
        assert sp.k1 == 0.0
        assert abs(sp.k2) == pytest.approx(1.0 / radius, rel=5e-3)
        assert abs(sp.d2 @ np.array([0.0, 0.0, 1.0])) < 1e-6


def test_a_planar_web_is_flat_with_exact_zeros():
    web = PartitionPatch(_cylinder(), "shear_web_0", np.linspace(0.0, 4.0, 5))
    sp = web.at(0.5, 1.3)
    assert sp.kind == "FLAT"
    assert (sp.k1, sp.k2, sp.theta) == (0.0, 0.0, 0.0)


def test_analytic_derivatives_match_finite_differences(cylinder_skin):
    patch = PartitionPatch(_twisted(31), SKIN_LE, np.linspace(0.0, 6.0, 31))
    t, z, eps = 0.4, 2.7, 1e-5
    d = patch.derivatives(t, z)
    assert np.allclose(d["x_t"], (patch(t + eps, z) - patch(t - eps, z)) / (2 * eps), atol=1e-5)
    assert np.allclose(d["x_z"], (patch(t, z + eps) - patch(t, z - eps)) / (2 * eps), atol=1e-5)
    xt = lambda zz: patch.derivatives(t, zz)["x_t"]  # noqa: E731
    assert np.allclose(d["x_tz"], (xt(z + eps) - xt(z - eps)) / (2 * eps), atol=1e-3)


def test_a_twisted_tapered_skin_is_doubly_curved():
    patch = PartitionPatch(_twisted(31), SKIN_LE, np.linspace(0.0, 6.0, 31))
    assert patch.sampling_ok
    assert patch.at(0.5, 3.0).kind == "DOUBLY"


def test_newton_projection_returns_the_foot_and_signed_distance(cylinder_skin):
    sp = cylinder_skin.at(0.45, 1.7)
    for offset in (0.01, -0.02):
        t, z, dist = cylinder_skin.closest(sp.x + offset * sp.normal, t0=0.3, z0=1.0)
        assert (t, z) == pytest.approx((0.45, 1.7), abs=1e-7)
        assert dist == pytest.approx(offset, abs=1e-9)


def test_sparse_stations_are_refused_for_a_curving_span():
    geo = _twisted(31)
    with pytest.raises(UnderSampledSpanLaw):
        PartitionPatch(geo, SKIN_LE, np.linspace(0.0, 6.0, 4))
    flagged = PartitionPatch(geo, SKIN_LE, np.linspace(0.0, 6.0, 4), sampling="flag")
    assert not flagged.sampling_ok


def test_bad_inputs():
    geo = _cylinder()
    with pytest.raises(ValueError, match="increasing"):
        PartitionPatch(geo, SKIN_LE, [1.0])
    with pytest.raises(ValueError, match="not on the shell graph"):
        PartitionPatch(geo, "no_such_run", [0.0, 1.0])
