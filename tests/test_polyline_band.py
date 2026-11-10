"""Numba unsigned polyline band matches the Python ``sdf_segment`` loop."""

from __future__ import annotations

import numpy as np

from implicit_geometry_kernel.sdf.kernels.polygon_numba import sdf_polyline_unsigned_numba
from implicit_geometry_kernel.sdf.kernels.primitives import sdf_segment
from implicit_geometry_kernel.sdf.nodes.tagged import PolylineBandSDF


def _python_polyline_unsigned(xy: np.ndarray, verts: np.ndarray) -> np.ndarray:
    x, y = xy[:, 0], xy[:, 1]
    d = np.full(xy.shape[0], np.inf)
    for a, b in zip(verts[:-1], verts[1:]):
        d = np.minimum(d, sdf_segment(x, y, a[0], a[1], b[0], b[1]))
    return np.asarray(d, dtype=np.float64)


def test_unsigned_polyline_numba_matches_segment_loop():
    verts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 0.5], [0.2, 0.8]], dtype=float)
    xx, yy = np.meshgrid(np.linspace(-0.2, 1.3, 25), np.linspace(-0.2, 1.1, 21))
    xy = np.column_stack([xx.ravel(), yy.ravel()])
    np.testing.assert_allclose(
        sdf_polyline_unsigned_numba(xy[:, 0], xy[:, 1], verts),
        _python_polyline_unsigned(xy, verts),
        atol=1e-12,
        rtol=1e-12,
    )


def test_polyline_band_inside_outside_radius():
    verts = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float)
    radius = 0.05
    band = PolylineBandSDF(verts, radius=radius, meta=3)
    pts = np.array([
        [0.5, 0.0, 0.0],
        [0.5, 0.02, 0.0],
        [0.5, 0.05, 0.0],
        [0.5, 0.2, 0.0],
    ], dtype=float)
    phi = band.evaluate(pts).phi
    assert phi[0] < 0.0
    assert phi[1] < 0.0
    assert abs(phi[2]) < 1e-12
    assert phi[3] > 0.0
    np.testing.assert_allclose(phi, _python_polyline_unsigned(pts[:, :2], verts) - radius, atol=1e-12)
