"""transform_second_moments_S_to_B must agree with integrating the rotated outline."""

import numpy as np
import pytest
from shapely.geometry import Polygon

from implicit_geometry_kernel import Airfoil
from implicit_geometry_kernel.frames import naca_to_S, rotate_S_to_B
from implicit_geometry_kernel.frames.tensors import (
    transform_second_moments_B_to_S,
    transform_second_moments_S_to_B,
)
from implicit_geometry_kernel.integration import section_properties_S


def _outline_S():
    a = Airfoil.naca4("2418")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 201)))
    v = np.vstack([np.column_stack([x[::-1], a.surface("upper", x)[::-1]]),
                   np.column_stack([x[1:], a.surface("lower", x)[1:]])])[:-1] * 2.0
    return naca_to_S(v, 2.0, 0.3)


@pytest.mark.parametrize("twist_deg", [-30.0, 7.5, 15.0, 45.0, 90.0])
def test_matches_rotated_outline(twist_deg):
    S = _outline_S()
    _, _, ixx, iyy, ixy = section_properties_S(Polygon(S))
    _, _, bxx, byy, bxy = section_properties_S(Polygon(rotate_S_to_B(S, np.radians(twist_deg))))
    got = transform_second_moments_S_to_B(ixx, iyy, ixy, twist_deg=twist_deg)
    np.testing.assert_allclose(got, (bxx, byy, bxy), rtol=1e-10, atol=1e-12)
    assert got[0] + got[1] == pytest.approx(ixx + iyy, rel=1e-12)
    back = transform_second_moments_B_to_S(*got, twist_deg=twist_deg)
    np.testing.assert_allclose(back, (ixx, iyy, ixy), rtol=1e-10, atol=1e-12)
