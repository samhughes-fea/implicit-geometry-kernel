"""Names and geometry that downstream code written against 0.2.x relies on."""

import pytest

import implicit_geometry_kernel as k
from implicit_geometry_kernel import Airfoil


def test_shear_web_names():
    from implicit_geometry_kernel.blade import SpanWeb
    from implicit_geometry_kernel.section.walls.shear_web import ShearWebComponent

    assert SpanWeb is k.ShearWeb
    assert k.ShearWeb is not ShearWebComponent


@pytest.mark.parametrize("code", ["0015", "0021"])
def test_naca4_n_is_points_per_surface(code):
    # ``n`` is points per surface, as in 0.2.0. At half that density the spline through
    # the points overshoots the nose and the leading edge opens by about 0.05% chord.
    a = Airfoil.naca4(code)
    gap = a.surface("upper", 0.0) - a.surface("lower", 0.0)
    assert abs(gap) < 2.5e-4
