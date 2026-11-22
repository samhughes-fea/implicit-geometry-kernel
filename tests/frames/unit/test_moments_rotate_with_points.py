"""CF-004: moment vectors rotate with the points that carry the load.

Axial forces F_z at points x_S give in-plane moments M_X = sum(Y F_z),
M_Y = -sum(X F_z) (R1.15). Rotating the points to B with ``rotate_S_to_B`` and
taking the moments there must equal ``transform_moments_S_to_B`` of the S-frame
moments (R1.16), and ``transform_moments_B_to_S`` must invert it. Before the fix
the moment transform rotated by the opposite twist.
"""

import numpy as np
import pytest

from implicit_geometry_kernel.frames import CoordinateFrameTransforms, rotate_S_to_B


def _moments(points, fz):
    return float((points[:, 1] * fz).sum()), float((-points[:, 0] * fz).sum())


@pytest.mark.parametrize("twist_deg", [-90.0, -35.0, -7.5, 0.0, 12.0, 20.0, 60.0, 90.0])
def test_moments_follow_rotated_points(twist_deg):
    rng = np.random.default_rng(4)
    pts_S = rng.normal(size=(50, 2))
    fz = rng.normal(size=50)
    mx_s, my_s = _moments(pts_S, fz)
    mx_b, my_b = _moments(rotate_S_to_B(pts_S, np.radians(twist_deg)), fz)
    t = CoordinateFrameTransforms()
    my_b_got, mx_b_got = t.transform_moments_S_to_B(my_s, mx_s, twist_deg)
    scale = abs(mx_s) + abs(my_s)
    assert abs(mx_b_got - mx_b) < 1e-12 * scale
    assert abs(my_b_got - my_b) < 1e-12 * scale
    my_s_back, mx_s_back = t.transform_moments_B_to_S(my_b, mx_b, twist_deg)
    assert abs(mx_s_back - mx_s) < 1e-12 * scale
    assert abs(my_s_back - my_s) < 1e-12 * scale


def test_moments_and_forces_rotate_like_points():
    t = CoordinateFrameTransforms()
    twist = 20.0
    pt_B = rotate_S_to_B(np.array([[1.0, 0.0]]), np.radians(twist))[0]
    vx, vy = t.transform_forces_S_to_B(1.0, 0.0, twist)
    my, mx = t.transform_moments_S_to_B(0.0, 1.0, twist)
    assert np.allclose([vx, vy], pt_B, atol=1e-15)
    assert np.allclose([mx, my], pt_B, atol=1e-15)
