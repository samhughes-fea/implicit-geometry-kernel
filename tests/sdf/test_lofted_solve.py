"""Lofted polygon SDF (exact at each z) and in-plane Newton projection."""

import numpy as np

from implicit_geometry_kernel.sdf import GeometrySolver, LoftedPolygonSDF
from implicit_geometry_kernel.sdf.nodes import CircleSDF, PolygonSDF


def _square(z):
    s = 1.0 + 0.5 * float(z)
    return np.array([[-s, -s], [s, -s], [s, s], [-s, s]], dtype=float)


def test_lofted_polygon_uses_own_station():
    field = LoftedPolygonSDF(_square, meta=4)
    inside_z0 = field.evaluate(np.array([[0.0, 0.0, 0.0]]))
    outside_z0 = field.evaluate(np.array([[1.2, 0.0, 0.0]]))
    inside_z2 = field.evaluate(np.array([[1.2, 0.0, 2.0]]))
    assert inside_z0.phi[0] < 0.0 and inside_z0.meta[0] == 4
    assert outside_z0.phi[0] > 0.0
    assert inside_z2.phi[0] < 0.0


def test_lofted_matches_fixed_polygon_at_same_z():
    z = 0.5
    lofted = LoftedPolygonSDF(_square)
    fixed = PolygonSDF(_square(z))
    rng = np.random.default_rng(0)
    xy = rng.uniform(-2.0, 2.0, size=(40, 2))
    pts = np.column_stack([xy, np.full(40, z)])
    np.testing.assert_allclose(lofted.evaluate(pts).phi, fixed.evaluate(pts).phi, atol=1e-12)


def test_project_circle_in_plane():
    solver = GeometrySolver(tolerance=1e-8, max_step=0.5)
    field = CircleSDF(0.0, 0.0, 1.0)
    start = np.array([[2.0, 0.0, 3.0], [0.0, 0.5, 3.0]])
    out = solver.project(field, start, plane=True)
    assert np.all(out.converged)
    np.testing.assert_allclose(np.linalg.norm(out.points[:, :2], axis=1), 1.0, atol=1e-6)
    np.testing.assert_allclose(out.points[:, 2], 3.0)


def test_intersect_two_circles():
    solver = GeometrySolver(tolerance=1e-8, max_step=0.5)
    a = CircleSDF(0.0, 0.0, 1.0)
    b = CircleSDF(1.0, 0.0, 1.0)
    guess = np.array([[0.5, 0.8, 0.0]])
    out = solver.intersect(a, b, guess)
    assert out.converged[0]
    np.testing.assert_allclose(out.points[0, 0], 0.5, atol=1e-5)
    np.testing.assert_allclose(abs(out.points[0, 1]), np.sqrt(0.75), atol=1e-5)
    assert out.points[0, 2] == 0.0
