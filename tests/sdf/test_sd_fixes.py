"""Regression tests for SD-001 to SD-004 (docs/sdf/bug_log.md, R4)."""

import numpy as np
import pytest

from implicit_geometry_kernel.extract import extract_zero_loops
from implicit_geometry_kernel.sdf import Box, Circle, Intersect, Offset, Union, eval_expr, sdf_ellipse
from implicit_geometry_kernel.sdf.kernels.polygon_bvh import build_polygon_accel


def _sampled_distance(points, rx, ry, n=400_000):
    t = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    boundary = np.column_stack([rx * np.cos(t), ry * np.sin(t)])
    d = np.empty(len(points))
    for k, p in enumerate(points):
        d[k] = np.min(np.hypot(*(boundary - p).T))
    inside = (points[:, 0] / rx) ** 2 + (points[:, 1] / ry) ** 2 <= 1.0
    return np.where(inside, -d, d)


@pytest.mark.parametrize("rx, ry", [(1.0, 0.5), (0.5, 1.0), (2.0, 0.3), (0.7, 0.7)])
def test_sd001_ellipse_is_the_distance(rx, ry):
    rng = np.random.default_rng(1)
    pts = np.vstack([rng.uniform(-2.0, 2.0, (60, 2)), [[0.0, 0.0], [0.3 * rx, 0.0], [0.0, 0.2 * ry],
                                                       [2.0 * rx, 0.0], [0.0, 2.0 * ry], [rx, 0.0]]])
    phi = sdf_ellipse(pts[:, 0], pts[:, 1], 0.0, 0.0, rx, ry)
    np.testing.assert_allclose(phi, _sampled_distance(pts, rx, ry), atol=1e-6)


def test_sd001_ellipse_closed_forms():
    assert sdf_ellipse(2.0, 0.0, rx=1.0, ry=0.5) == pytest.approx(1.0, abs=1e-14)
    assert sdf_ellipse(0.0, 0.0, rx=2.0, ry=1.0) == pytest.approx(-1.0, abs=1e-14)
    assert sdf_ellipse(1.5, 0.25, cx=1.5, cy=-0.25, rx=1.0, ry=1.0) == pytest.approx(-0.5, abs=1e-12)


def test_sd002_offset_of_disjoint_intersection_box_contains_region():
    e = Offset(Intersect(Box(0, 0, 1, 1), Box(2.5, 0, 1, 1)), 0.5)
    xs, ys = np.meshgrid(np.linspace(-3, 6, 451), np.linspace(-3, 3, 301))
    inside = eval_expr(e, xs, ys) <= 0.0
    assert inside.any()
    x0, x1, y0, y1 = e.bbox()
    assert x0 <= xs[inside].min() and xs[inside].max() <= x1
    assert y0 <= ys[inside].min() and ys[inside].max() <= y1


def test_sd003_two_close_circles_extract_as_two_loops():
    two = Union(Circle(-0.37, -0.37, 0.5), Circle(0.37, 0.37, 0.5))
    loops = extract_zero_loops(two, 0.0, (-1.5, 1.5, -1.5, 1.5), n=48)
    assert len(loops) == 2
    for loop in loops:
        np.testing.assert_allclose(loop[0], loop[-1], atol=1e-9)
        assert np.abs(two.evaluate(loop).phi).max() < 1e-7


def test_sd004_accel_holds_edges_and_box_only():
    square = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    acc = build_polygon_accel(square)
    assert set(acc) == {"ax", "ay", "bx", "by", "aabb"}
    np.testing.assert_array_equal(acc["bx"], [1.0, 1.0, 0.0, 0.0])
    assert acc["aabb"] == (0.0, 1.0, 0.0, 1.0)
