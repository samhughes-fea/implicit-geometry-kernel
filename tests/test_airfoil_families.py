"""CST, t/c families, flatback TE, Selig/Lednicer readers, interpolate by t/c."""

from pathlib import Path

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil


def test_naca5_and_naca6():
    five = Airfoil.naca5("23012")
    assert five.max_thickness == pytest.approx(0.12, rel=0.08)
    assert float(np.max(five.surface("upper", np.linspace(0.0, 1.0, 80)))) > 0.02
    six = Airfoil.naca6(63, 15)
    assert six.max_thickness == pytest.approx(0.15, rel=0.08)


def test_cst_symmetric_closed():
    foil = Airfoil.cst([0.2, 0.2, 0.2], [0.2, 0.2, 0.2], te_thickness=0.0)
    x = np.linspace(0.0, 1.0, 50)
    np.testing.assert_allclose(foil.surface("upper", x), -foil.surface("lower", x), atol=1e-9)
    assert foil.max_thickness > 0.05


def test_thickness_family_and_flatback():
    base = Airfoil.naca4("0012")
    thick = base.with_thickness(0.18)
    assert thick.max_thickness == pytest.approx(0.18, rel=1e-3)
    x = np.linspace(0.0, 1.0, 40)
    camber = 0.5 * (base.surface("upper", x) + base.surface("lower", x))
    camber2 = 0.5 * (thick.surface("upper", x) + thick.surface("lower", x))
    np.testing.assert_allclose(camber[1:], camber2[1:], atol=1e-7)
    blunt = base.with_te_thickness(0.02)
    te = float(blunt.surface("upper", 1.0) - blunt.surface("lower", 1.0))
    assert te == pytest.approx(0.02, abs=1e-6)


def test_interpolate_by_thickness():
    a, b = Airfoil.naca4("0012"), Airfoil.naca4("0024")
    mid = Airfoil.interpolate_by_thickness(a, b, 0.18)
    assert mid.max_thickness == pytest.approx(0.18, rel=1e-3)


def test_selig_and_lednicer_readers(tmp_path: Path):
    foil = Airfoil.naca4("2412", n=41)
    selig = foil.contour_points(41)
    selig_path = tmp_path / "naca2412.dat"
    selig_path.write_text("NACA 2412\n" + "\n".join(f"{x:.8f} {y:.8f}" for x, y in selig), encoding="utf-8")
    back = Airfoil.from_selig(selig_path)
    x = np.linspace(0.05, 0.9, 20)
    np.testing.assert_allclose(back.surface("upper", x), foil.surface("upper", x), atol=5e-3)

    x_le = np.linspace(0.0, 1.0, 21)
    yu, yl = foil.surface("upper", x_le), foil.surface("lower", x_le)
    led_path = tmp_path / "lednicer.dat"
    lines = ["NACA 2412", "21. 21."]
    lines += [f"{x:.8f} {y:.8f}" for x, y in zip(x_le, yu)]
    lines += [f"{x:.8f} {y:.8f}" for x, y in zip(x_le, yl)]
    led_path.write_text("\n".join(lines), encoding="utf-8")
    led = Airfoil.from_lednicer(led_path)
    np.testing.assert_allclose(led.surface("upper", x), foil.surface("upper", x), atol=5e-3)
