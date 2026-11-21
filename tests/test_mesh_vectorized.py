"""Vectorized smooth, triangle packing, and quality match the pointwise path."""

import numpy as np

from implicit_geometry_kernel.mesh.fill import (
    _edge_clearance,
    _smooth_interior,
    triangulate_ring,
)
from implicit_geometry_kernel.mesh.fill import points_in_ring
from implicit_geometry_kernel.mesh.quality import (
    QUALITY_DTYPE,
    _cell_area,
    _corner_angles,
    _element_metrics,
    _equiangle_skew,
    _volume_ratio,
    _warp,
    element_quality,
)
from implicit_geometry_kernel.mesh.shell import _weld_triangles


def _smooth_pointwise(cloud, tris, n_boundary, xy, h):
    """The per-vertex neighbour mean, kept as the oracle for the array smooth."""
    if n_boundary >= len(cloud):
        return None
    nbrs = [[] for _ in range(len(cloud))]
    for a, b, c in tris:
        for u, v in ((int(a), int(b)), (int(b), int(c)), (int(c), int(a))):
            nbrs[u].append(v)
            nbrs[v].append(u)
    moved = np.asarray(cloud[n_boundary:], float).copy()
    changed = False
    margin = 0.35 * float(h)
    for i in range(n_boundary, len(cloud)):
        nb = nbrs[i]
        if len(nb) < 2:
            continue
        acc = cloud[np.unique(nb)].mean(axis=0)
        if not points_in_ring(xy, acc.reshape(1, 2))[0]:
            continue
        if float(_edge_clearance(xy, acc.reshape(1, 2))[0]) < margin:
            continue
        moved[i - n_boundary] = acc
        changed = True
    if not changed:
        return None
    return moved


def _cloud():
    boundary = np.array([[0.0, 0.0], [1.0, 0.0], [1.2, 0.8], [0.1, 1.0]], float)
    interior = np.array([[0.35, 0.30], [0.70, 0.45], [0.40, 0.70], [0.85, 0.25]], float)
    cloud = np.vstack([boundary, interior])
    tris = np.array([
        [0, 1, 4], [1, 7, 4], [1, 2, 7], [2, 5, 7],
        [2, 3, 5], [3, 6, 5], [3, 0, 6], [0, 4, 6],
        [4, 7, 5], [4, 5, 6],
    ], np.int32)
    return cloud, tris, len(boundary), boundary, 0.35


def test_smooth_matches_the_pointwise_mean():
    cloud, tris, n_boundary, xy, h = _cloud()
    new = _smooth_interior(cloud, tris, n_boundary, xy, h)
    old = _smooth_pointwise(cloud, tris, n_boundary, xy, h)
    assert new is not None and old is not None
    assert np.allclose(new, old)


def test_triangulate_ring_is_one_triangle_array():
    ring = np.array([
        [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    ])
    faces = triangulate_ring(ring, size=0.35, pattern="unstructured", split_boundary=False)
    listed = [faces[i] for i in range(len(faces))]
    assert faces.shape == (len(listed), 3, 3)
    assert np.allclose(faces, np.stack(listed))
    assert len(faces) >= 2


def _quality_rows(nodes, elements, thickness):
    """The per-element quality loop, kept as the oracle for the array path."""
    nodes = np.asarray(nodes, float)
    thickness = np.asarray(thickness, float)
    table = np.zeros(len(elements), dtype=QUALITY_DTYPE)
    area = np.zeros(len(elements), float)
    for i, element in enumerate(elements):
        p = nodes[list(element)]
        angles = _corner_angles(p)
        aspect, jac = _element_metrics(p)
        cell_area = _cell_area(p)
        area[i] = cell_area
        table[i] = (
            aspect,
            _equiangle_skew(angles),
            float(np.min(angles)) if len(angles) else 0.0,
            float(np.max(angles)) if len(angles) else 0.0,
            float(np.min(jac)) if len(jac) else 0.0,
            cell_area,
            0.0,
            1.0,
            _warp(p),
        )
    volume = area * thickness
    table["shell_volume"] = volume
    table["volume_ratio"] = _volume_ratio(elements, volume)
    return table


def test_quality_rows_match_the_element_loop():
    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.2, 0.8, 0.0],
        [1.0, 0.0, 0.0],
        [1.1, 0.9, 0.05],
        [0.0, 1.0, 0.0],
        [-0.1, 0.2, 0.0],
        [2.0, 0.0, 0.1],
        [2.2, 1.0, 0.0],
        [1.5, 1.2, 0.3],
        [1.4, 0.2, -0.1],
    ], float)
    elements = [
        (0, 1, 2),
        (3, 4, 5, 6),
        (7, 8, 9, 10),
        (2, 5, 6),
        (0, 1, 5, 2),
    ]
    thickness = [0.02, 0.04, 0.01, 0.03, 0.05]
    got = element_quality(nodes, elements, thickness)
    expect = _quality_rows(nodes, elements, thickness)
    for name in QUALITY_DTYPE.names:
        assert np.allclose(got[name], expect[name], equal_nan=True)


def test_station_weld_appends_new_points_in_unique_order():
    nodes = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
    origins = ["station", "station"]
    index = {(0.0, 0.0, 0.0): 0, (1.0, 0.0, 0.0): 1}
    faces = np.array([
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.2, 0.3, 0.0]],
        [[1.0, 0.0, 0.0], [0.5, 0.4, 0.0], [0.2, 0.3, 0.0]],
        [[0.2, 0.3, 0.0], [0.2, 0.3, 0.0], [0.5, 0.4, 0.0]],
    ])
    corners = _weld_triangles(faces, index, nodes, origins)
    assert corners.shape == (3, 3)
    assert list(corners[0]) == [0, 1, 3]
    assert list(corners[1]) == [1, 2, 3]
    assert list(corners[2]) == [-1, -1, -1]
    assert origins[2:] == ["interior", "interior"]
    assert [tuple(p) for p in nodes[2:]] == [(0.5, 0.4, 0.0), (0.2, 0.3, 0.0)]
