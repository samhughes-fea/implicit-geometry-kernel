"""Native Steiner Delaunay fill of a planar n-gon."""

import numpy as np

import pytest

from implicit_geometry_kernel.mesh.fill import triangulate_ring


def test_triangulate_square_is_triangles():
    ring = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        float,
    )
    faces = triangulate_ring(ring, size=0.25)
    assert len(faces)
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
    assert len(faces)
    assert all(len(f) == 3 for f in faces)
    for face in faces:
        c = face[:, :2].mean(axis=0)
        in_bay = (c[0] < 2.0 - 1e-6) and (1.0 + 1e-6 < c[1] < 2.0 - 1e-6)
        assert not in_bay


def _edge_of(faces, a, b, tol=1e-8) -> bool:
    for face in faces:
        for i in range(3):
            p, q = face[i], face[(i + 1) % 3]
            if (np.linalg.norm(p - a) < tol and np.linalg.norm(q - b) < tol) or (
                np.linalg.norm(p - b) < tol and np.linalg.norm(q - a) < tol
            ):
                return True
    return False


def _on_segment(p, a, b, tol=1e-6) -> bool:
    ab = b - a
    lab = float(np.dot(ab, ab))
    if lab < 1e-24:
        return float(np.linalg.norm(p - a)) <= tol
    t = float(np.dot(p - a, ab) / lab)
    if t < -1e-6 or t > 1.0 + 1e-6:
        return False
    return float(np.linalg.norm(p - (a + t * ab))) <= tol


def _segment_covered(faces, a, b, tol=1e-6) -> bool:
    """Triangle edges that lie on ``a``–``b`` cover the whole segment."""
    a, b = np.asarray(a, float)[:2], np.asarray(b, float)[:2]
    ab = b - a
    lab = float(np.dot(ab, ab))
    intervals = []
    for face in faces:
        for i in range(3):
            p, q = face[i, :2], face[(i + 1) % 3, :2]
            if not (_on_segment(p, a, b, tol) and _on_segment(q, a, b, tol)):
                continue
            t0 = float(np.dot(p - a, ab) / lab)
            t1 = float(np.dot(q - a, ab) / lab)
            lo, hi = (t0, t1) if t0 <= t1 else (t1, t0)
            if hi - lo > 1e-8:
                intervals.append((lo, hi))
    if not intervals:
        return False
    intervals.sort()
    cursor = 0.0
    for lo, hi in intervals:
        if lo > cursor + 1e-4:
            return False
        cursor = max(cursor, hi)
    return cursor >= 1.0 - 1e-4


def test_unstructured_uses_ring_edges_and_one_size():
    # Four segments a side, already the target length, so none are split.
    side = np.linspace(0.0, 1.0, 5)
    ring = []
    ring.extend([x, 0.0, 0.0] for x in side[:-1])
    ring.extend([1.0, y, 0.0] for y in side[:-1])
    ring.extend([x, 1.0, 0.0] for x in side[::-1][:-1])
    ring.extend([0.0, y, 0.0] for y in side[::-1][:-1])
    ring = np.asarray(ring, float)
    faces = triangulate_ring(ring, size=0.25, pattern="unstructured")
    assert len(faces)
    for i, a in enumerate(ring):
        assert _edge_of(faces, a, ring[(i + 1) % len(ring)])
    lengths = []
    for face in faces:
        for i in range(3):
            lengths.append(float(np.linalg.norm(face[(i + 1) % 3] - face[i])))
    lengths = np.asarray(lengths)
    assert lengths.min() > 0.08
    assert lengths.max() < 0.7
    covered = sum(
        abs(float(np.dot(f[:, 0], np.roll(f[:, 1], -1)) - np.dot(f[:, 1], np.roll(f[:, 0], -1)))) / 2.0
        for f in faces
    )
    assert abs(covered - 1.0) < 0.05


def test_unstructured_splits_long_edges():
    ring = np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    faces = triangulate_ring(ring, size=0.3, pattern="unstructured")
    for i, a in enumerate(ring):
        b = ring[(i + 1) % len(ring)]
        assert _segment_covered(faces, a, b)
        assert not _edge_of(faces, a, b)
    lengths = [
        float(np.linalg.norm(face[(i + 1) % 3] - face[i]))
        for face in faces for i in range(3)
    ]
    assert max(lengths) < 0.75
    covered = sum(
        abs(float(np.dot(f[:, 0], np.roll(f[:, 1], -1)) - np.dot(f[:, 1], np.roll(f[:, 0], -1)))) / 2.0
        for f in faces
    )
    assert abs(covered - 4.0) < 0.15


def test_unstructured_can_keep_the_ring_vertices_as_the_boundary():
    ring = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    faces = triangulate_ring(ring, size=0.3, pattern="unstructured", split_boundary=False)
    verts = np.unique(np.round(np.vstack(faces), 8), axis=0)

    def on_boundary(p) -> bool:
        return min(p[0], 2.0 - p[0], p[1], 1.0 - p[1]) < 1e-6

    boundary = [p for p in verts if on_boundary(p)]
    assert len(boundary) == 4
    for p in boundary:
        assert any(np.linalg.norm(p - a) < 1e-6 for a in ring)
    for i, a in enumerate(ring):
        assert _segment_covered(faces, a, ring[(i + 1) % len(ring)])


def test_unstructured_keeps_the_notch():
    ring = np.array([
        [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 3.0, 0.0], [0.0, 3.0, 0.0],
        [0.0, 2.0, 0.0], [2.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    ])
    faces = triangulate_ring(ring, size=0.45, pattern="unstructured")
    for i, a in enumerate(ring):
        assert _segment_covered(faces, a, ring[(i + 1) % len(ring)])
    for face in faces:
        c = face[:, :2].mean(axis=0)
        in_bay = (c[0] < 2.0 - 1e-6) and (1.0 + 1e-6 < c[1] < 2.0 - 1e-6)
        assert not in_bay
