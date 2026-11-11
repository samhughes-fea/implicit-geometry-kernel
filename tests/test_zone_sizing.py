"""Zone sizing: node counts follow region length, constant within a zone, odd, still conforming."""

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
from implicit_geometry_kernel.blade.topology.zones import _odd_nodes, zone_sized_templates, zones_of


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


def test_zone_sized_strips_are_quads_and_share_station_nodes():
    geo = _geo()
    cuts, tmpls, smps = [], [], []
    for z in np.linspace(0.0, 10.0, 11):
        cut = section_cut_near(geo, float(z))
        t = _template_from_cut(geo.region_graph, cut, float(z), 8)
        cuts.append(cut)
        tmpls.append(t)
        smps.append(sample_region_polylines(cut, cut.shell_graph(), t, geo.region_graph))
    h = 0.05
    sized = zone_sized_templates(tmpls, smps, h)
    smps = [sample_region_polylines(c, c.shell_graph(), t, geo.region_graph) for c, t in zip(cuts, sized)]
    for zone in zones_of(sized):
        for name in smps[zone[0]]:
            counts = {len(smps[i][name]) for i in zone if name in smps[i]}
            assert len(counts) == 1 and counts.pop() % 2 == 1, name
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
        zone_sized_templates([], [], 0.0)
