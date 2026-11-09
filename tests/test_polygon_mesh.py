"""Kernel polygonal fill of a planar n-gon (no Gmsh)."""

import numpy as np

from implicit_geometry_kernel.mesh.polygon import point_in_ring, points_in_ring, polygon_mesh_ring


def test_points_in_ring_matches_scalar():
    ring = np.array([[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0],
                     [0.0, 2.0], [2.0, 2.0], [2.0, 1.0], [0.0, 1.0]])
    pts = np.array([[1.5, 0.5], [1.0, 1.5], [2.5, 2.5], [0.0, 0.0]])
    batched = points_in_ring(ring, pts)
    scalar = np.array([point_in_ring(p, ring) for p in pts])
    assert np.array_equal(batched, scalar)
    assert batched[0] and not batched[1] and batched[2]


def test_polygon_mesh_square_has_small_ngons():
    ring = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        float,
    )
    faces = polygon_mesh_ring(ring, size=0.25, lloyd=0)
    assert len(faces) >= 4
    parent = ring[:, :2]
    areas = []
    for face in faces:
        assert len(face) > 4
        xy = face[:, :2]
        areas.append(abs(float(np.dot(xy[:, 0], np.roll(xy[:, 1], -1)) - np.dot(xy[:, 1], np.roll(xy[:, 0], -1)))) / 2.0)
        assert point_in_ring(xy.mean(axis=0), parent)
    assert max(areas) < 0.35


def test_polygon_mesh_square_covers_boundary():
    ring = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        float,
    )
    faces = polygon_mesh_ring(ring, size=0.25, lloyd=0)
    edges = []
    for face in faces:
        xy = face[:, :2]
        for i in range(len(xy)):
            edges.append((xy[i], xy[(i + 1) % len(xy)]))
    covered = 0
    parent = ring[:, :2]
    for i in range(4):
        a, b = parent[i], parent[(i + 1) % 4]
        ab = b - a
        hit = False
        for p, q in edges:
            mid = 0.5 * (p + q)
            t = np.dot(mid - a, ab) / (np.dot(ab, ab) + 1e-30)
            proj = a + np.clip(t, 0.0, 1.0) * ab
            if float(np.linalg.norm(mid - proj)) < 2e-3 and 0.02 < t < 0.98:
                hit = True
                break
        covered += int(hit)
    assert covered == 4


def test_polygon_mesh_stays_inside_c_cell():
    ring = np.array([
        [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 3.0, 0.0], [0.0, 3.0, 0.0],
        [0.0, 2.0, 0.0], [2.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    ])
    parent = ring[:, :2]
    faces = polygon_mesh_ring(ring, size=0.45, lloyd=0)
    assert faces
    for face in faces:
        assert len(face) > 4
        c = face[:, :2].mean(axis=0)
        assert point_in_ring(c, parent)
        in_bay = (c[0] < 2.0 - 1e-6) and (1.0 + 1e-6 < c[1] < 2.0 - 1e-6)
        assert not in_bay
