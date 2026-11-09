"""Native Steiner Delaunay fill of a planar n-gon."""

import numpy as np

import pytest

from implicit_geometry_kernel.mesh.delaunay import triangulate_ring


def test_triangulate_square_is_triangles():
    ring = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        float,
    )
    faces = triangulate_ring(ring, size=0.25)
    assert faces
    assert all(len(f) == 3 for f in faces)
    for face in faces:
        xy = face[:, :2]
        area = abs(float(np.dot(xy[:, 0], np.roll(xy[:, 1], -1)) - np.dot(xy[:, 1], np.roll(xy[:, 0], -1)))) / 2.0
        assert area < 0.15
        c = xy.mean(axis=0)
        assert 0.0 - 1e-6 <= c[0] <= 1.0 + 1e-6
        assert 0.0 - 1e-6 <= c[1] <= 1.0 + 1e-6
    assert abs(sum(
        abs(float(np.dot(f[:, 0], np.roll(f[:, 1], -1)) - np.dot(f[:, 1], np.roll(f[:, 0], -1)))) / 2.0
        for f in faces
    ) - 1.0) < 0.08


@pytest.mark.parametrize("engine", [None, "numba"])
def test_triangulate_stays_inside_c_cell(engine):
    ring = np.array([
        [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 3.0, 0.0], [0.0, 3.0, 0.0],
        [0.0, 2.0, 0.0], [2.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    ])
    faces = triangulate_ring(ring, size=0.45, engine=engine)
    assert faces
    assert all(len(f) == 3 for f in faces)
    for face in faces:
        c = face[:, :2].mean(axis=0)
        in_bay = (c[0] < 2.0 - 1e-6) and (1.0 + 1e-6 < c[1] < 2.0 - 1e-6)
        assert not in_bay
