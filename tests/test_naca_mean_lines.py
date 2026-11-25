"""NACA 5-digit and 6-series mean lines against their definitions, and the morph TE gap (BG-055 to BG-057)."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil
from implicit_geometry_kernel.blade.loft.loft_tables import _morph_surfaces, _morph_surfaces_many, _smoothstep
from implicit_geometry_kernel.sdf.shapes.naca import (
    naca_five_digit_vertices,
    six_series_camber_thickness,
)

# Abbott & von Doenhoff, normal 5-digit mean lines at C_l,i = 0.3: P -> (r, k1)
ABBOTT = {1: (0.0580, 361.400), 2: (0.1260, 51.640), 3: (0.2025, 15.957), 4: (0.2900, 6.643), 5: (0.3910, 3.230)}


def _pairs(verts):
    n = len(verts) // 2
    return verts[: n + 1][::-1], verts[n:]


@pytest.mark.parametrize("P", sorted(ABBOTT))
@pytest.mark.parametrize("L", [1, 2, 4])
def test_five_digit_mean_line_is_abbotts(P, L):
    r, k1 = ABBOTT[P]
    k1 *= (3.0 * L / 20.0) / 0.3
    up, lo = _pairs(naca_five_digit_vertices(L, 10 * P, 12, 201, 1.0))
    xc = 0.5 * (up[:, 0] + lo[:, 0])
    yc = np.where(xc <= r, k1 / 6 * (xc**3 - 3 * r * xc**2 + r * r * (3 - r) * xc), k1 * r**3 / 6 * (1 - xc))
    np.testing.assert_allclose(0.5 * (up[:, 1] + lo[:, 1]), yc, atol=1e-15)


def test_naca_23012_camber_and_position():
    up, lo = _pairs(naca_five_digit_vertices(2, 30, 12, 801, 1.0))
    mid = 0.5 * (up + lo)
    k = int(np.argmax(mid[:, 1]))
    assert mid[k, 1] == pytest.approx(0.01839, abs=2e-5)
    assert mid[k, 0] == pytest.approx(0.15, abs=5e-3)


def test_q_digit_selects_the_reflexed_mean_line():
    np.testing.assert_allclose(
        naca_five_digit_vertices(2, 31, 12, 121, 1.0),
        naca_five_digit_vertices(2, 30, 12, 121, 1.0, reflex=True),
    )
    up, lo = _pairs(naca_five_digit_vertices(2, 31, 12, 401, 1.0))
    assert abs(0.5 * (up[-1, 1] + lo[-1, 1])) < 1e-12  # reflexed line returns to the chord at the TE


@pytest.mark.parametrize("pq", [32, 60, 0])
def test_bad_middle_digits_raise(pq):
    with pytest.raises(ValueError):
        naca_five_digit_vertices(2, pq, 12, 121, 1.0)


@pytest.mark.parametrize("family", [63, 64])
def test_six_series_keeps_camber_and_scales_thickness(family):
    x = np.linspace(0.02, 0.98, 97)
    ref_camber, _ = six_series_camber_thickness(family, x)
    for t in (12, 18, 24):
        foil = Airfoil.naca6(family, t, n=201)
        np.testing.assert_allclose(0.5 * (foil.surface("upper", x) + foil.surface("lower", x)), ref_camber, atol=2e-4)
        assert foil.max_thickness == pytest.approx(t / 100.0, rel=1e-3)
    assert Airfoil.naca6(64, 18).name == "NACA 64-218"


def test_smoothstep_morph_moves_the_te_gap_with_the_blend_weight():
    circle, flat = Airfoil.ellipse(1.0), Airfoil.naca4("2418").with_te_thickness(0.03)
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 161)))
    ws = np.array([0.25, 0.5, 0.75])
    up_many, lo_many = _morph_surfaces_many(circle, flat, ws, x)
    for i, w in enumerate(ws):
        up, lo = _morph_surfaces(circle, flat, float(w), x)
        assert up[-1] - lo[-1] == pytest.approx(_smoothstep(float(w)) * 0.03, abs=1e-12)
        np.testing.assert_allclose(up, up_many[i], atol=1e-14)
        np.testing.assert_allclose(lo, lo_many[i], atol=1e-14)
