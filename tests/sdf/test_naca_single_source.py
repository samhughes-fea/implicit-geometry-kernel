"""NACA vertex maths lives in ``sdf.shapes.naca``; Airfoil is the shape, AirfoilSDF the field."""

from __future__ import annotations

import inspect

from implicit_geometry_kernel.blade.airfoil import Airfoil
from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF
from implicit_geometry_kernel.sdf.shapes.naca import naca_four_digit_vertices


def test_naca4_and_from_naca_call_the_shared_generator() -> None:
    assert "naca_four_digit_vertices" in inspect.getsource(Airfoil.naca4)
    assert "naca_four_digit_vertices" in inspect.getsource(AirfoilSDF.from_naca)
    verts = naca_four_digit_vertices(0.02, 0.4, 0.12, 81, 1.0)
    assert verts.ndim == 2 and verts.shape[1] == 2


def test_to_sdf_returns_airfoil_sdf() -> None:
    field = Airfoil.naca4("2412").to_sdf(chord=2.0, n=81)
    assert isinstance(field, AirfoilSDF)
    assert field.chord == 2.0
