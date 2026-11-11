"""Region counts: node counts follow region length (element size h/2 to h), odd, still conforming."""

from collections import Counter

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.topology import (
    _template_from_cut,
    loft_mapped_strips,
    sample_region_polylines,
    section_cut_near,
)
from implicit_geometry_kernel.blade.topology.map_loft import _transition_quads
from implicit_geometry_kernel.blade.topology.region_counts import _held_edges, _odd_nodes, region_counts


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.004, pivot_fraction=0.3,
                         spar_cap_thickness=0.008, shear_web_thickness=0.004)


def test_odd_nodes_follow_length():
    assert _odd_nodes(1.0, 0.25) == 5
    assert _odd_nodes(1.1, 0.25) == 7
    assert _odd_nodes(0.01, 0.25) == 3
    assert all(_odd_nodes(x, 0.1) % 2 == 1 for x in np.linspace(0.0, 2.0, 37))


def test_counts_are_held_while_the_element_size_stays_within_h_over_2_and_h():
    h = 0.1
    assert _held_edges(1.0, None, h, 3) == 10
    assert _held_edges(0.9, 10, h, 3) == 10  # 0.09 per element: held
    assert _held_edges(0.5, 10, h, 3) == 10  # h / 2 exactly: held
    assert _held_edges(0.49, 10, h, 3) == 6  # below h / 2: from length
    assert _held_edges(1.01, 10, h, 3) == 12  # past h: from length
    assert _held_edges(0.0, 4, h, 3) == 2


def test_sized_strips_are_quads_and_share_station_nodes():
    geo = _geo()
    cuts, tmpls, smps = [], [], []
    for z in np.linspace(0.0, 10.0, 11):
        cut = section_cut_near(geo, float(z))
        t = _template_from_cut(geo.region_graph, cut, float(z), 8)
        cuts.append(cut)
        tmpls.append(t)
        smps.append(sample_region_polylines(cut, cut.shell_graph(), t, geo.region_graph))
    h = 0.05
    sized = region_counts(tmpls, smps, h)
    smps = [sample_region_polylines(c, c.shell_graph(), t, geo.region_graph) for c, t in zip(cuts, sized)]
    for smp in smps:
        for name, pts in smp.items():
            length = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
            edges = len(pts) - 1
            assert edges % 2 == 0, name
            assert length / edges <= h * (1 + 1e-9), name
            assert edges == 2 or length / edges >= 0.5 * h * (1 - 1e-9), name
    # counts follow length: the longest region gets more nodes than the shortest
    lengths = {n: float(np.sum(np.linalg.norm(np.diff(p, axis=0), axis=1))) for n, p in smps[0].items()}
    assert len(smps[0][max(lengths, key=lengths.get)]) > len(smps[0][min(lengths, key=lengths.get)])
    edges: Counter = Counter()
    for t0, s0, t1, s1 in zip(sized, smps, sized[1:], smps[1:]):
        for _name, face, _c in loft_mapped_strips(s0, s1, t0, t1):
            assert len(face) == 4
            ks = [tuple(np.round(p, 9)) for p in face]
            for a, b in zip(ks, ks[1:] + ks[:1]):
                if a != b:
                    edges[tuple(sorted((a, b)))] += 1
    for (a, b), n in edges.items():
        if n == 1:
            assert a[2] == b[2] and a[2] in (0.0, 10.0)


def test_zone_sizing_rejects_bad_input():
    with pytest.raises(ValueError):
        region_counts([], [], 0.0)


def test_transition_zip_is_the_same_whichever_station_is_finer():
    def worst(na, nb):
        pa = np.column_stack([np.linspace(0, 1, na), np.zeros(na), np.zeros(na)])
        pb = np.column_stack([np.linspace(0, 1, nb), np.zeros(nb), np.full(nb, 2.0 / (max(na, nb) - 1))])
        out = []
        for _n, f, _c in _transition_quads("x", pa, pb):
            for i in range(4):
                a, b = f[i - 1] - f[i], f[(i + 1) % 4] - f[i]
                out.append(np.degrees(np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b))))
        return min(out), max(out)

    for na, nb in ((3, 5), (5, 7), (7, 9)):
        assert worst(na, nb) == pytest.approx(worst(nb, na))
        lo, hi = worst(nb, na)
        assert lo >= 45.0 - 1e-6 and hi <= 143.2
