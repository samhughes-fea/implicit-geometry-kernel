"""The polygon fill conforms to its ring: no hanging node against the skin."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.topology import (
    _template_from_cut,
    rib_cell_polygons_from_samples,
    sample_region_polylines,
)
from implicit_geometry_kernel.mesh.polygon import _area2, _ring_basis, polygon_mesh_ring


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
                         spar_cap_thickness=0.008, shear_web_thickness=0.004)


def _rings(geo, z, member, n_chord):
    cut = geo.section_at(z)
    tmpl = _template_from_cut(geo.region_graph, cut, z, n_chord)
    samples = sample_region_polylines(cut, cut.shell_graph(), tmpl, geo.region_graph)
    return rib_cell_polygons_from_samples(tmpl, samples, member)


def _check(ring):
    ring = np.asarray(ring, float)
    faces = polygon_mesh_ring(ring)
    assert faces
    origin, u, v, xy = _ring_basis(ring)
    area = 0.5 * abs(_area2(xy))
    ring_nodes = {tuple(p) for p in ring}
    a, ab = xy, np.roll(xy, -1, axis=0) - xy
    total = 0.0
    for face in faces:
        assert len(face) > 4
        fxy = np.column_stack([(face - origin) @ u, (face - origin) @ v])
        assert _area2(fxy) > 0.0
        total += 0.5 * _area2(fxy)
        for p, q in zip(face, fxy):
            if tuple(p) in ring_nodes:
                continue
            t = np.clip(((q - a) * ab).sum(axis=1) / (ab * ab).sum(axis=1), 0.0, 1.0)
            # a node that is not a ring node is off the ring: no hanging node
            assert np.min(np.linalg.norm(a + t[:, None] * ab - q, axis=1)) > 1e-9 * np.sqrt(area)
    assert {tuple(p) for f in faces for p in f} >= ring_nodes
    assert total == pytest.approx(area, rel=1e-12)


@pytest.mark.parametrize("n_chord", [8, 16])
def test_tip_closure_fill_conforms_to_the_skin(n_chord):
    rings = _rings(_geo(), 10.0, "tip_closure", n_chord)
    assert len(rings) == 3
    for _name, ring in rings:
        _check(ring)


def test_rib_fill_conforms_to_the_skin():
    geo = _geo()
    rib = geo.transverse[1]
    for _name, ring in _rings(geo, rib.z, rib, 8):
        _check(ring)


def test_square_with_coarse_ring_still_conforms():
    # four ring nodes cannot carry cells a quarter of the side: the fill coarsens
    ring = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    faces = polygon_mesh_ring(ring, size=0.1)
    assert faces
    for face in faces:
        for p in face:
            on_edge = np.isclose(p[0], 0.0) or np.isclose(p[0], 1.0) or np.isclose(p[1], 0.0) or np.isclose(p[1], 1.0)
            assert not on_edge or any(np.array_equal(p, r) for r in ring)
