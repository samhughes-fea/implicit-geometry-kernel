"""Batched mesh tables, area vectors, and loop projections match the pointwise path."""

import numpy as np

from implicit_geometry_kernel.blade.topology.map_loft import _span_mids
from implicit_geometry_kernel.mesh.quality import _corner_angles
from implicit_geometry_kernel.mesh.shell import (
    _area_vector,
    _areas_means,
    _closest_on_loop,
    _closest_on_loop_many,
    _oriented,
)
from implicit_geometry_kernel.mesh.tables import (
    _corner_angles_block,
    _diagonal_min,
    _element_normals,
    _face_normal,
    diagonal_mins,
)


def _faces():
    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.2, 0.0],
        [0.2, 1.0, 0.1],
        [1.1, 0.9, 0.4],
        [0.0, 0.0, 1.0],
        [2.0, 0.0, 0.2],
        [2.2, 1.0, 0.0],
        [1.4, 1.3, 0.5],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.2, 0.0],
    ], float)
    elements = [
        (0, 1, 2),
        (1, 3, 4, 2),
        (4, 5, 6, 7),
        (0, 1, 3, 7, 2),
        (8, 9),
    ]
    return nodes, elements


def test_normals_angles_and_diagonals_match_the_element_loop():
    nodes, elements = _faces()
    normals = _element_normals(nodes, elements)
    expect = [
        _face_normal(nodes[list(element)]) if len(element) >= 3 else np.zeros(3)
        for element in elements
    ]
    assert np.allclose(normals, expect)

    for n in (3, 4):
        chosen = [nodes[list(element)] for element in elements if len(element) == n]
        block = _corner_angles_block(np.stack(chosen))
        for row, face in zip(block, chosen):
            assert np.allclose(row, _corner_angles(face))

    class _Mesh:
        pass

    mesh = _Mesh()
    mesh.nodes = nodes
    mesh.elements = elements
    got = diagonal_mins(mesh)
    expect_d = tuple(_diagonal_min(nodes, element) for element in elements)
    assert np.allclose(got, expect_d)


def test_area_vectors_match_the_polygon_normal():
    nodes, elements = _faces()
    indices = [i for i, element in enumerate(elements) if len(element) >= 3]
    areas, means = _areas_means(nodes, elements, indices)
    for row, mean, i in zip(areas, means, indices):
        face = nodes[list(elements[i])]
        assert np.allclose(row, _area_vector(face))
        assert np.allclose(mean, face.mean(axis=0))


def test_skin_orientation_survives_a_centroid_nearer_another_station():
    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 1.0],
        [0.0, 0.0, 1.0],
    ], float)
    out, flipped = _oriented(nodes, [(0, 1, 2, 3)], ["skin"], [0.0, 0.4, 1.0], {})
    assert out == [(0, 1, 2, 3)]
    assert flipped == [0]


def test_inward_skin_sheet_flips_outward():
    nodes = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [0.0, 0.0, 1.0],
        [1.0, 1.0, 0.0], [1.0, 1.0, 1.0],
        [0.0, 1.0, 0.0], [0.0, 1.0, 1.0],
    ], float)
    # Each face is wound so its area vector points into the square tube.
    elements = [
        (0, 3, 2, 1),
        (1, 2, 5, 4),
        (4, 5, 7, 6),
        (6, 7, 3, 0),
    ]
    out, flipped = _oriented(nodes, elements, ["skin"] * 4, [0.0, 1.0], {})
    assert flipped == [1, 1, 1, 1]
    for element in out:
        area = _area_vector(nodes[list(element)])
        centre = nodes[list(element)].mean(axis=0) - np.array([0.5, 0.5, 0.5])
        assert float(area[:2] @ centre[:2]) > 0.0


def _closest_pointwise(loop, pts):
    return np.stack([_closest_on_loop(loop, p) for p in pts])


def test_closest_points_match_one_query_at_a_time():
    loop = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    pts = np.array([
        [0.2, -0.4, 0.0],
        [1.0, 0.3, 0.0],
        [2.4, 0.5, 0.0],
        [0.5, 1.2, 0.1],
    ])
    assert np.allclose(_closest_on_loop_many(loop, pts), _closest_pointwise(loop, pts))


def test_span_mids_project_the_interior_together():
    pa = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    pb = pa + np.array([0.0, 0.0, 2.0])
    calls = []

    def on_section(_name, xyz, at_end):
        calls.append((np.asarray(xyz, float).shape, bool(at_end)))
        xyz = np.asarray(xyz, float)
        if xyz.ndim == 1:
            return xyz + np.array([0.0, 1.0, 0.0])
        return xyz + np.array([0.0, 2.0, 0.0])

    mids = _span_mids("skin", pa, pb, on_section)
    assert calls[0] == ((3,), True)
    assert calls[1] == ((3,), True)
    assert calls[2] == ((2, 3), False)
    assert np.allclose(mids[0, 1], 1.0)
    assert np.allclose(mids[-1, 1], 1.0)
    assert np.allclose(mids[1:-1, 1], 2.0)


def test_off_station_seam_blends_the_station_keys(monkeypatch):
    import importlib

    mesh_mod = importlib.import_module("implicit_geometry_kernel.mesh.shell.shell_mesh")

    def refuse(*_args, **_kwargs):
        raise AssertionError("off-station seam built a section")

    monkeypatch.setattr(mesh_mod, "section_cut_near", refuse)

    def square(z, shift):
        return np.array([
            [0.0, shift, z],
            [1.0, shift, z],
            [1.0, 1.0 + shift, z],
            [0.0, 1.0 + shift, z],
        ], float)

    class Cut:
        def __init__(self, z, shift):
            self.z = z
            self.chord, self.twist, self.pitch_axis, self.offset = 1.0, 0.0, 0.0, (0.0, 0.0)
            self.loop = square(z, shift)
            self.keys = {
                ("le", "suction"): np.array([0.0, shift, z]),
                ("te", "suction"): np.array([5.0, 5.0, z]),
            }

        def shell_graph(self):
            return self

    left, right = Cut(0.0, 0.0), Cut(10.0, 2.0)

    def skin(_graph, z):
        return _graph.loop

    def keys(cut, _graph, _region, _z):
        return cut.keys

    monkeypatch.setattr(mesh_mod, "_skin_mid_loop", skin)
    monkeypatch.setattr(mesh_mod, "place_station_keys", keys)

    class Geo:
        region_graph = object()
        _mesh_z_cuts = {0.0: left, 10.0: right}

        @staticmethod
        def interpolate(_z):
            return 1.0, 0.0, 0.0, (0.0, 0.0)

        @staticmethod
        def _placement_twist(twist):
            return float(twist)

    project, _loop_at = mesh_mod._section_projector(Geo(), set())
    snapped = project("", np.array([0.1, 1.0, 5.0]), True)
    assert np.allclose(snapped, (0.0, 1.0, 5.0))
