"""Mapped loft strips meet at the station nodes: no cracks along a station line."""

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


def _geo(key):
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, key, 0.004, pivot_fraction=0.3,
                         spar_cap_thickness=0.008, shear_web_thickness=0.004)


def _key(p):
    return tuple(np.round(np.asarray(p, float), 9))


@pytest.mark.parametrize("key", ["0A-LC-0R", "2D-F-LC-0R"])
def test_strips_share_the_station_nodes(key):
    geo = _geo(key)
    samples = []
    for z in np.linspace(0.0, 10.0, 6):
        cut = section_cut_near(geo, float(z))
        tmpl = _template_from_cut(geo.region_graph, cut, float(z), 8)
        samples.append((float(z), tmpl, sample_region_polylines(cut, cut.shell_graph(), tmpl, geo.region_graph)))
    edges: Counter = Counter()
    for (_za, t0, s0), (_zb, t1, s1) in zip(samples, samples[1:]):
        for _name, face, _collapsing in loft_mapped_strips(s0, s1, t0, t1):
            ks = [_key(p) for p in face]
            for a, b in zip(ks, ks[1:] + ks[:1]):
                if a != b:
                    edges[tuple(sorted((a, b)))] += 1
    for z, _t, smp in samples[1:-1]:
        station = {_key(p) for pts in smp.values() for p in np.asarray(pts, float)}
        on_line = {p for e in edges for p in e if abs(p[2] - z) < 1e-9}
        # every strip node on an interior station line is one of that station's nodes
        assert on_line == station
    # strips are open only along the root and tip lines (the tip closure closes the tip)
    free = [e for e, n in edges.items() if n == 1]
    assert free
    for a, b in free:
        assert a[2] == b[2] and a[2] in (0.0, 10.0)
