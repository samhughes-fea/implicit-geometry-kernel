"""Cap keys sit on the cut's own cap edges, and a cap bonds across a concave skin vertex."""

import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from implicit_geometry_kernel.section.junctions.cap_skin_bonds import (
    _Arc,
    _bond_one_cap,
    _spread_shared_feet,
)

ROOT = Path(__file__).resolve().parents[1]


def test_shared_feet_are_spread_by_cap_arc_length():
    sigma = np.array([0.0, 1.0, 1.0, 1.0, 3.0])
    prefix = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    out = _spread_shared_feet(sigma, prefix, "cap")
    assert np.all(np.diff(out) > 0.0)
    assert out[0] == 0.0 and out[-1] == 3.0
    assert out[1:4] == pytest.approx([0.75, 1.5, 2.25])


def test_shared_feet_at_a_cap_end_are_refused():
    with pytest.raises(ValueError):
        _spread_shared_feet(np.array([1.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]), "cap")


def test_cap_bonds_round_a_concave_skin_vertex():
    # the skin bends towards the cap's far side: two cap vertices on the offset
    # arc round the corner both have the corner as their nearest skin point
    skin = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 1.0]])
    r = 0.1
    cap = np.array([[0.2, -r], [0.6, -r], [1.0, -r],
                    [1.0 + r * np.sin(np.pi / 4), -r * np.cos(np.pi / 4)],
                    [1.6 + r * np.sin(np.pi / 4), 0.6 - r * np.cos(np.pi / 4)]])
    strip = SimpleNamespace(label="skin", midline_twisted_pitch_axis_blade_frame=skin)
    chain = [_Arc(0, strip, reversed=False)]
    new_cap, _inserts, bonds = _bond_one_cap(
        SimpleNamespace(label="cap", midline_twisted_pitch_axis_blade_frame=cap), chain)
    assert len(bonds) == 1
    host = np.asarray(bonds[0].host_points)
    # one-to-one: the skin partners advance strictly along the skin
    seg = np.diff(host, axis=0)
    assert np.all(np.linalg.norm(seg, axis=1) > 0.0)
    assert np.all(seg @ np.array([1.0, 0.5]) > 0.0)
    # no geometry moves: every skin partner lies on the skin polyline
    for p in host:
        a, b = skin[:-1], skin[1:]
        t = np.clip(((p - a) * (b - a)).sum(1) / ((b - a) ** 2).sum(1), 0.0, 1.0)
        assert np.min(np.linalg.norm(a + t[:, None] * (b - a) - p, axis=1)) < 1e-12
    assert len(new_cap) >= len(cap)


def _catalogue():
    env = os.environ.get("REFERENCE_BLADES_DIR")
    roots = ([Path(env)] if env else []) + [ROOT / "reference_blades", ROOT.parent / "reference-blades"]
    for path in roots:
        if (path / "wind" / "iea22" / "IEA-22-280-RWT.yaml").is_file() or (
            path / "iea22" / "IEA-22-280-RWT.yaml"
        ).is_file():
            return path
    return None


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_iea22_cap_keys_are_the_cut_cap_edges():
    import runpy

    from implicit_geometry_kernel.blade.topology import section_cut_near
    from implicit_geometry_kernel.blade.topology.map import place_station_keys
    from implicit_geometry_kernel.section.assembly.members import wall_type_for_label

    ns = runpy.run_path(str(ROOT / "examples" / "blade" / "08_reference_blades" / "reference_blades.py"))
    _name, geo, _ = ns["load_windio"](ns["catalogue_yaml"](_catalogue(), "iea22/IEA-22-280-RWT.yaml"))

    def is_cap(label):
        wt = wall_type_for_label(label)
        return wt is not None and wt.kind == "cap"

    for z in (6.89, 12.4, 13.78, 15.2, 27.56):
        cut = section_cut_near(geo, z)
        graph = cut.shell_graph()
        # the pressure cap bonds to its skin near the root (a concave skin vertex)
        assert {"spar_cap_suction_surface", "spar_cap_pressure_surface"} <= {r.name for r in graph.loop}, z
        edges = [np.asarray(e.point, float)[:2] for e in graph.partition_edges()
                 if is_cap(e.name_a) != is_cap(e.name_b)]
        keys = place_station_keys(cut, graph, geo.region_graph, float(cut.z))
        cap_keys = [p for (name, _side), p in keys.items() if str(name).startswith("cap_")]
        assert cap_keys, z
        for p in cap_keys:
            assert min(np.linalg.norm(e - p[:2]) for e in edges) < 1e-9, z


def _iea34_catalogue():
    cat = _catalogue()
    if cat is None:
        return None
    for rel in ("wind/iea3.4_iea10/IEA-3.4-130-RWT.yaml", "iea3.4_iea10/IEA-3.4-130-RWT.yaml"):
        if (cat / rel).is_file():
            return cat
    return None


@pytest.mark.skipif(_iea34_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_iea34_cap_edge_passing_a_web_gets_a_station_and_stays_conforming():
    import runpy
    from collections import Counter

    from implicit_geometry_kernel.blade.topology import (
        _template_from_cut,
        loft_mapped_strips,
        sample_region_polylines,
        section_cut_near,
    )
    from implicit_geometry_kernel.blade.topology.crossings import key_crossing_stations

    ns = runpy.run_path(str(ROOT / "examples" / "blade" / "08_reference_blades" / "reference_blades.py"))
    _name, geo, _ = ns["load_windio"](ns["catalogue_yaml"](_iea34_catalogue(), "iea3.4_iea10/IEA-3.4-130-RWT.yaml"))
    zs = [3.15, 4.725]
    extra = key_crossing_stations(geo, zs)
    assert len(extra) == 1 and zs[0] < extra[0] < zs[1]
    stations = []
    for z in sorted(zs + list(extra)):
        cut = section_cut_near(geo, z)
        tmpl = _template_from_cut(geo.region_graph, cut, z, 8)
        stations.append((tmpl, sample_region_polylines(cut, cut.shell_graph(), tmpl, geo.region_graph)))
    webs = {w.name for t, _s in stations for w in t.webs}
    edges: Counter = Counter()
    for (t0, s0), (t1, s1) in zip(stations, stations[1:]):
        for name, face, _c in loft_mapped_strips(s0, s1, t0, t1):
            if name in webs:
                continue
            assert len(face) == 4, name
            ks = [tuple(np.round(p, 9)) for p in face]
            for a, b in zip(ks, ks[1:] + ks[:1]):
                if a != b:
                    edges[tuple(sorted((a, b)))] += 1
    # skin strips are open only along the first and last station lines
    for (a, b), n in edges.items():
        if n == 1:
            assert a[2] == b[2] and a[2] in (zs[0], zs[1])
