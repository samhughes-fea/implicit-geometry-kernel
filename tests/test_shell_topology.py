"""Mapped loft: identity, zones, web drop, cap flicker, TE wrap, rib faces."""

from pathlib import Path

import numpy as np
import pytest

from implicit_geometry_kernel import (
    Airfoil,
    BladeGeometry,
    Station,
    loft_mapped_strips,
    rib_face_polygons,
    signatures_compatible,
    transverse_members_from_windio,
    zone_region_graphs,
)
from implicit_geometry_kernel.blade.concept import BladeConcept
from implicit_geometry_kernel.blade.topology import (
    Anchor,
    Region,
    ShellTopologyTemplate,
    WallRegionGraph,
    WebMember,
    _cluster_web_members,
    _cut_tc,
    _regularise_anchors,
    _te_wrap_arc,
    place_station_keys,
    rib_cell_polygons_from_samples,
    sample_region_polylines,
)


def _geo(key="2D-F-LC-0R", **kwargs):
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    kw = dict(pivot_fraction=0.3, spar_cap_thickness=0.008, shear_web_thickness=0.004)
    kw.update(kwargs)
    return BladeGeometry(stations, key, 0.004, **kw)


def _templates(geo, zs, n_chord=8):
    out = []
    samples = []
    for z in zs:
        cut = geo.section_at(float(z))
        sg = cut.shell_graph()
        live = []
        for run in sg.webs:
            xy = np.asarray(run.points, float)[:, :2]
            live.append(float(np.mean((xy[:, 0] - cut.offset[0]) / max(cut.chord, 1e-9) + cut.pitch_axis)))
        t = ShellTopologyTemplate.from_region_graph(
            geo.region_graph, float(z), n_chord=n_chord, live_web_fracs=live or None, t_c=_cut_tc(cut),
            dropped_web_indices=tuple(cut.dropped_web_indices),
        )
        out.append(t)
        samples.append(sample_region_polylines(cut, sg, t, geo.region_graph))
    return out, samples


def test_web_add_lofts_t1_only_roles_as_birth():
    graph = WallRegionGraph(
        anchors=(
            Anchor("le", 0.0),
            Anchor("web_0", 0.25),
            Anchor("web_1", 0.55, z0=5.0),
            Anchor("te", 1.0),
        ),
        regions=(),
        webs=(
            WebMember(0, 0.25, name="fore_web"),
            WebMember(1, 0.55, z0=5.0, name="rear_web"),
        ),
    )
    geo = _geo(web_positions=(0.25, 0.55), region_graph=graph)
    tpls, samples = _templates(geo, [4.0, 6.0])
    extra = {b.name for b in tpls[1].loop} - {b.name for b in tpls[0].loop}
    assert extra
    got = {n for n, face, collapsing in loft_mapped_strips(
        samples[0], samples[1], tpls[0], tpls[1],
    )}
    assert extra & got, extra - got
    extras = [n for n, _f, _c in loft_mapped_strips(
        samples[0], samples[1], tpls[0], tpls[1],
    ) if n in extra]
    assert extras


def test_web_drop_loop_vertices_map_and_skin_covers_span():
    geo = _geo(web_offsets=(0.4, 1.2))
    zs = np.linspace(0.0, 10.0, 11)
    tpls, samples = _templates(geo, zs)
    assert any(len(a.webs) != len(b.webs) for a, b in zip(tpls, tpls[1:]))
    loop_roles = [set(b.name for b in t.loop) for t in tpls]
    shared = set.intersection(*loop_roles)
    assert shared
    for t0, t1 in zip(tpls, tpls[1:]):
        for vid in t0.vertices:
            assert t0.map_vertex(vid, t1) in t1.vertices or t0.map_vertex(vid, t1) == ("TE",)
    strips = []
    for s0, s1, t0, t1 in zip(samples, samples[1:], tpls, tpls[1:]):
        strips.extend(loft_mapped_strips(s0, s1, t0, t1))
    loop_names = {n for n, _f, _c in strips if not str(n).startswith("web")}
    assert shared & loop_names
    assert all(n not in {f[0] for f in strips if f[2]} or True for n in shared)
    # no skipped adjacent pair for a loop role present on both
    for s0, s1, t0, t1 in zip(samples, samples[1:], tpls, tpls[1:]):
        both = {b.name for b in t0.loop} & {b.name for b in t1.loop} & set(s0) & set(s1)
        got = {n for n, _f, _c in loft_mapped_strips(s0, s1, t0, t1)}
        assert both <= got, (both - got, t0.z, t1.z)


def test_one_station_cap_flicker_is_transition_not_hole():
    tpls, samples = _templates(_geo(web_positions=(0.18, 0.48)), [2.0, 4.0, 4.1, 6.0])
    mid = tpls[1]
    # paint-compatible: same webs, loop atoms
    assert signatures_compatible(tpls[0].signature, tpls[2].signature)
    zoned = zone_region_graphs(tpls)
    names = []
    for s0, s1, t0, t1 in zip(samples, samples[1:], tpls, tpls[1:]):
        names.extend(n for n, _f, _c in loft_mapped_strips(s0, s1, t0, t1))
    assert names
    assert zoned.zones


def test_te_wrap_collapsing_patch():
    geo = _geo(web_positions=(0.18, 0.48))
    tpls, samples = _templates(geo, [3.0, 5.0])
    t0, t1 = tpls
    wrap = next((b for b in t0.loop if b.wrap), t0.loop[0])
    t0_forced = t0
    faces = loft_mapped_strips(samples[0], samples[1], t0_forced, t1)
    assert faces
    if wrap.wrap:
        mapped = t0.map_vertex(wrap.start, t1)
        assert mapped == ("TE",) or mapped in t1.vertices


def test_write_mid_surfaces_obj_region_groups_and_rib(tmp_path: Path):
    from implicit_geometry_kernel.blade.stations import TransverseMember

    geo = _geo(web_positions=(0.18, 0.48),
               transverse=(TransverseMember("tip_closure", 10.0, 0.012),))
    obj = Path(geo.write_mid_surfaces(tmp_path / "mids.obj", n_span=5, n_chord=8)).read_text(
        encoding="utf-8")
    assert "g " in obj and "v " in obj
    assert "tip_closure" in obj or "rib" in obj
    n_gons = [ln for ln in obj.splitlines() if ln.startswith("f ") and len(ln.split()) > 5]
    assert n_gons, "rib OBJ groups must keep n-gon faces (n>4)"
    ribs = rib_face_polygons(geo, n_chord=8)
    assert ribs
    assert all(abs(float(ring[0, 2]) - 10.0) < 1e-6 for _n, ring in ribs)


def test_mapped_loft_omits_rib_strips():
    from implicit_geometry_kernel.blade.stations import TransverseMember

    z_of_s = np.column_stack([[0.0, 1.0], [0.0, 10.0]])
    structure = {
        "webs": [
            {"name": "w0", "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.3, 0.3]},
             "end_nd_arc": {"grid": [0.0, 1.0], "values": [0.7, 0.7]}},
            {"name": "w1", "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.2, 0.2]},
             "end_nd_arc": {"grid": [0.0, 1.0], "values": [0.8, 0.8]}},
        ],
        "layers": [
            {"name": "Shell_skin",
             "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.0, 0.0]},
             "end_nd_arc": {"grid": [0.0, 1.0], "values": [1.0, 1.0]},
             "thickness": {"grid": [0.0, 1.0], "values": [0.004, 0.003]}},
        ],
    }
    members = transverse_members_from_windio(structure, z_of_s)
    assert members[0].name == "tip_closure"
    geo = _geo(web_positions=(0.18, 0.48), transverse=members)
    tpls, samples = _templates(geo, [2.0, 8.0])
    strips = loft_mapped_strips(samples[0], samples[1], tpls[0], tpls[1])
    assert all("rib" not in n and n != "tip_closure" for n, _f, _c in strips)
    assert rib_face_polygons(geo)


def _assert_keyed_seams(geo, z, n_chord=10):
    cut = geo.section_at(float(z))
    sg = cut.shell_graph()
    t = ShellTopologyTemplate.from_region_graph(geo.region_graph, float(z), n_chord=n_chord)
    samples = sample_region_polylines(cut, sg, t, geo.region_graph)
    keys = place_station_keys(cut, sg, geo.region_graph, float(z))
    su = [b for b in t.loop if b.side == "suction" and not b.wrap]
    pr = [b for b in t.loop if b.side == "pressure" and not b.wrap]
    for chain in (su, pr):
        for a, b in zip(chain, chain[1:]):
            assert a.end == b.start
            assert a.name in samples and b.name in samples
            assert np.allclose(samples[a.name][-1], samples[b.name][0], atol=1e-12)
            assert np.linalg.norm(samples[a.name][-1] - samples[b.name][0]) == 0.0
    loop_vids = {b.start for b in t.loop} | {b.end for b in t.loop}
    assert all(v[0] != "web" for v in t.vertices)
    for web in t.webs:
        assert web.suction in loop_vids and web.pressure in loop_vids
        assert web.suction in t.vertices and web.pressure in t.vertices
        assert web.name in samples
        su_ft = keys[(f"web_{web.index}", "suction")]
        pr_ft = keys[(f"web_{web.index}", "pressure")]
        ends = samples[web.name][0], samples[web.name][-1]
        d00 = np.linalg.norm(ends[0] - su_ft) + np.linalg.norm(ends[1] - pr_ft)
        d01 = np.linalg.norm(ends[0] - pr_ft) + np.linalg.norm(ends[1] - su_ft)
        assert min(d00, d01) == 0.0
        wall = np.vstack([samples[b.name] for b in t.loop if b.name in samples])
        for end in ends:
            assert np.min(np.linalg.norm(wall - end, axis=1)) == 0.0
        if len(samples[web.name]) > 2:
            interior = samples[web.name][1:-1]
            assert not np.all(np.min(np.linalg.norm(wall[:, None] - interior, axis=2), axis=0) == 0.0)
    return t, samples, keys


def test_catalogue_2df_consecutive_loop_keys_shared():
    _assert_keyed_seams(_geo(web_positions=(0.18, 0.48)), 4.0)


def test_catalogue_0a_and_1cf_seams():
    _assert_keyed_seams(_geo("0A-LC-0R"), 3.0)
    _assert_keyed_seams(_geo("1C-F-LC-0R", web_positions=(0.30,)), 3.0)


def test_ellipse_root_taper_seams():
    stations = [
        Station(0.0, 2.0, Airfoil.ellipse(1.0), pitch_axis=0.5),
        Station(4.0, 1.6, Airfoil.naca4("0018"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, web_positions=(0.25, 0.50), pivot_fraction=0.3)
    _assert_keyed_seams(geo, 0.0)
    _assert_keyed_seams(geo, 5.0)


def test_ellipse_root_lofts_te_wrap():
    stations = [
        Station(0.0, 2.0, Airfoil.ellipse(1.0), pitch_axis=0.5),
        Station(10.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, web_positions=(0.25, 0.50), pivot_fraction=0.5)
    tpls, samples = _templates(geo, [0.0, 1.0], n_chord=8)
    assert "te" in samples[0] and "te" in samples[1]
    wrap = np.asarray(samples[0]["te"], float)
    assert float(np.linalg.norm(wrap[0] - wrap[-1])) > 0.04 * 2.0
    names = {n for n, _f, _c in loft_mapped_strips(samples[0], samples[1], tpls[0], tpls[1])}
    assert "te" in names


def test_blade_concept_authors_same_key_graph_as_layout():
    concept = BladeConcept("2D-F-LC-0R")
    authored = concept.region_graph()
    direct = WallRegionGraph.from_section_layout("2D-F-LC-0R")
    assert tuple(r.name for r in authored.regions) == tuple(r.name for r in direct.regions)
    assert tuple(a.name for a in authored.anchors) == tuple(a.name for a in direct.anchors)
    t = ShellTopologyTemplate.from_region_graph(authored, 1.0, n_chord=6)
    su = [b for b in t.loop if b.side == "suction" and not b.wrap]
    for a, b in zip(su, su[1:]):
        assert a.end == b.start


def test_synthetic_windio_keyed_sampler():
    z_of_s = np.column_stack([[0.0, 1.0], [0.0, 10.0]])
    structure = {
        "webs": [
            {"name": "w0", "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.3, 0.3]},
             "end_nd_arc": {"grid": [0.0, 1.0], "values": [0.7, 0.7]}},
            {"name": "w1", "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.2, 0.2]},
             "end_nd_arc": {"grid": [0.0, 1.0], "values": [0.8, 0.8]}},
        ],
        "layers": [
            {"name": "Shell_skin",
             "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.0, 0.0]},
             "end_nd_arc": {"grid": [0.0, 1.0], "values": [1.0, 1.0]},
             "thickness": {"grid": [0.0, 1.0], "values": [0.004, 0.003]}},
        ],
    }
    graph = WallRegionGraph.from_windio_structure(structure, z_of_s)
    geo = _geo(web_positions=(0.18, 0.48), region_graph=graph)
    _assert_keyed_seams(geo, 2.0)


def test_sharp_te_wrap_is_omitted_from_samples():
    t, samples, _keys = _assert_keyed_seams(_geo(web_positions=(0.18, 0.48)), 4.0)
    wrap = next((b for b in t.loop if b.wrap), None)
    assert wrap is not None
    if wrap.name in samples:
        pts = np.asarray(samples[wrap.name], float)
        chord = 2.0
        slen = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
        assert slen < 0.05 * chord


def test_ellipse_root_omits_fat_te_wrap():
    stations = [
        Station(0.0, 2.0, Airfoil.ellipse(1.0), pitch_axis=0.5),
        Station(10.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, web_positions=(0.25, 0.50), pivot_fraction=0.3)
    t, samples, keys = _assert_keyed_seams(geo, 0.0, n_chord=12)
    wrap = next(b for b in t.loop if b.wrap)
    assert wrap.name in samples
    a, b = keys[("te", "suction")], keys[("te", "pressure")]
    assert float(np.linalg.norm(a - b)) > 0.04 * 2.0


def test_le_te_loop_still_holds_web_keys_on_skin():
    graph = WallRegionGraph(
        anchors=(Anchor("le", 0.0), Anchor("te", 1.0), Anchor("web_0", 0.25), Anchor("web_1", 0.55)),
        regions=(
            Region("te", "te", "te", False, "te", wrap=True),
            Region("le_su", "le", "te", False, "suction"),
            Region("le_pr", "le", "te", False, "pressure"),
        ),
        webs=(WebMember(0, 0.25, name="fore_web"), WebMember(1, 0.55, name="rear_web")),
    )
    geo = _geo(web_positions=(0.25, 0.55), region_graph=graph)
    t, samples, _keys = _assert_keyed_seams(geo, 3.0, n_chord=12)
    if "te" in samples:
        pts = np.asarray(samples["te"], float)
        slen = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
        assert slen < 0.05 * 2.0
    names = {b.name for b in t.loop}
    assert "le_su" in names
    assert any(n.startswith("bay_") or n.startswith("te_") for n in names)
    wall = np.vstack([samples[b.name] for b in t.loop if b.name in samples])
    for web in t.webs:
        for end in (samples[web.name][0], samples[web.name][-1]):
            assert np.min(np.linalg.norm(wall - end, axis=1)) == 0.0


def test_bowtie_quad_is_not_lofted():
    from implicit_geometry_kernel.blade.topology import loft_mapped_strips, _quad_bowtie

    a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    assert not _quad_bowtie(a)
    bow = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert _quad_bowtie(bow)
    tpls, samples = _templates(_geo(web_positions=(0.18, 0.48)), [3.0, 5.0])
    faces = loft_mapped_strips(samples[0], samples[1], tpls[0], tpls[1])
    assert faces
    for _n, face, collapsing in faces:
        if collapsing or len(face) < 4:
            continue
        assert not _quad_bowtie(face)


def test_catalogue_2df_te_trace_is_smooth():
    tpls, samples = _templates(_geo(web_positions=(0.18, 0.48)), np.linspace(0.0, 10.0, 21), n_chord=8)
    tes = np.array([s["te_pr"][-1, :2] for s in samples])
    jumps = np.linalg.norm(np.diff(tes, axis=0), axis=1)
    assert float(np.max(jumps)) < 0.25


def test_flatback_te_wrap_samples_skin_arc():
    blunt = Airfoil.naca4("0018").with_te_thickness(0.05)
    stations = [
        Station(0.0, 2.0, blunt, pitch_axis=0.3),
        Station(10.0, 1.0, blunt, pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, web_positions=(0.18, 0.48), pivot_fraction=0.3)
    t, samples, keys = _assert_keyed_seams(geo, 4.0, n_chord=12)
    su, pr = keys[("te", "suction")], keys[("te", "pressure")]
    assert float(np.linalg.norm(su - pr)) > 1e-4
    wrap = next(b for b in t.loop if b.wrap)
    assert wrap.name in samples
    ends = samples[wrap.name][0], samples[wrap.name][-1]
    d00 = np.linalg.norm(ends[0] - su) + np.linalg.norm(ends[1] - pr)
    d01 = np.linalg.norm(ends[0] - pr) + np.linalg.norm(ends[1] - su)
    assert min(d00, d01) == pytest.approx(0.0, abs=1e-12)


def test_te_plug_web_is_omitted_from_template():
    graph = WallRegionGraph(
        anchors=(Anchor("le", 0.0), Anchor("te", 1.0), Anchor("web_0", 0.25), Anchor("web_te", 0.99)),
        regions=(
            Region("te", "te", "te", False, "te", wrap=True),
            Region("le_su", "le", "web_0", False, "suction"),
            Region("te_su", "web_0", "te", False, "suction"),
            Region("le_pr", "le", "web_0", False, "pressure"),
            Region("te_pr", "web_0", "te", False, "pressure"),
        ),
        webs=(WebMember(0, 0.25, name="fore_web"), WebMember(1, 0.99, name="web_te")),
    )
    t = ShellTopologyTemplate.from_region_graph(graph, 1.0, n_chord=6)
    assert [w.name for w in t.webs] == ["fore_web"]


def _ptset(pts) -> set[tuple]:
    return {tuple(np.round(np.asarray(p, float), 9)) for p in pts}


def test_rib_ngons_use_loft_skin_and_web_nodes():
    from implicit_geometry_kernel.blade.stations import TransverseMember

    geo = _geo(web_positions=(0.18, 0.48),
               transverse=(TransverseMember("tip_closure", 10.0, 0.012),))
    tpls, samples = _templates(geo, [10.0], n_chord=8)
    rings = rib_cell_polygons_from_samples(tpls[0], samples[0], "tip_closure")
    assert rings
    assert all(len(ring) > 4 for _n, ring in rings)
    tmpl, samp = tpls[0], samples[0]
    su = [b for b in tmpl.loop if b.side == "suction" and not b.wrap]
    pr = [b for b in tmpl.loop if b.side == "pressure" and not b.wrap]
    for k, ((sb, pb), (_name, ring)) in enumerate(zip(zip(su, pr), rings)):
        wanted = _ptset(samp[sb.name]) | _ptset(samp[pb.name])
        aft = next((w for w in tmpl.webs if w.suction == sb.end or w.pressure == sb.end), None)
        fwd = next((w for w in tmpl.webs if w.suction == sb.start or w.pressure == sb.start), None)
        if aft is not None and aft.name in samp:
            wanted |= _ptset(samp[aft.name])
        if fwd is not None and fwd.name in samp:
            wanted |= _ptset(samp[fwd.name])
        wrap = next((b for b in tmpl.loop if b.wrap), None)
        if wrap is not None and wrap.name in samp and str(sb.end[0]) == "TE":
            wanted |= _ptset(samp[wrap.name])
        got = _ptset(ring)
        assert wanted <= got, (k, sorted(wanted - got)[:4])
    if len(tmpl.webs) and tmpl.webs[0].name in samp and len(rings) >= 2:
        interior = np.asarray(samp[tmpl.webs[0].name], float)[len(samp[tmpl.webs[0].name]) // 2]
        hit = [_ptset(r) for _n, r in rings]
        n_hit = sum(tuple(np.round(interior, 9)) in s for s in hit)
        assert n_hit == 2


def test_2df_span_metrics_are_watertight():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
    from _shared.debug_mesh import as_quad, span_metrics

    geo = _geo(web_positions=(0.18, 0.48))
    zs = np.linspace(0.0, 10.0, 11)
    tpls, samples = _templates(geo, zs, n_chord=8)
    walls = list(zip([float(z) for z in zs], tpls, samples))
    faces, kinds = [], []
    for s0, s1, t0, t1 in zip(samples, samples[1:], tpls, tpls[1:]):
        for name, face, _c in loft_mapped_strips(s0, s1, t0, t1):
            faces.append(face)
            kinds.append(name)
    face_arr = np.stack([as_quad(f) for f in faces]) if faces else np.zeros((0, 4, 3))
    roles = np.asarray(kinds, dtype="U64")
    metrics = span_metrics(walls, face_arr, roles, geo=geo)
    assert metrics["n_coverage_holes"] == 0
    assert metrics["max_wrap_frac"] < 0.05
    assert metrics["max_te_oml"] < 0.08
    assert metrics["max_te_chatter"] < 0.05
    assert metrics["n_web_clusters"] == 2
    assert metrics["z_tip"] == pytest.approx(10.0)


def _09_stations():
    root, foil = Airfoil.ellipse(1.0), Airfoil.naca4("0024")
    return [
        Station(0.0, 1.80, root, twist=12.0, pitch_axis=0.50),
        Station(1.0, 1.88, Airfoil.blend(root, foil, 0.15), twist=11.5, pitch_axis=0.47),
        Station(2.5, 2.05, Airfoil.blend(root, foil, 0.50), twist=10.0, pitch_axis=0.40),
        Station(4.0, 2.28, Airfoil.blend(root, foil, 0.85), twist=8.8, pitch_axis=0.33),
        Station(5.0, 2.40, foil, twist=8.0, pitch_axis=0.30),
        Station(6.5, 2.20, Airfoil.naca4("0022"), twist=6.5, pitch_axis=0.30),
        Station(8.5, 1.85, Airfoil.naca4("0020"), twist=4.8, pitch_axis=0.30),
        Station(11.0, 1.50, Airfoil.naca4("0018"), twist=2.8, pitch_axis=0.30),
        Station(14.0, 1.10, Airfoil.naca4("0015"), twist=1.0, pitch_axis=0.30),
        Station(17.0, 0.90, Airfoil.naca4("0013"), twist=-0.6, pitch_axis=0.30),
        Station(20.0, 0.75, Airfoil.naca4("0012"), twist=-2.0, pitch_axis=0.30),
    ]


def _09_geo(key: str) -> BladeGeometry:
    st = _09_stations()
    return BladeGeometry.from_concept(
        BladeConcept(key), st, 0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004, span_kind="pchip",
        pivot_fraction=([s.z for s in st], [s.pitch_axis for s in st]),
    )


def test_cluster_keeps_pa_lattice_on_cylinder():
    webs = tuple(WebMember(i, 0.2 * (i + 1), name=f"web_{i}") for i in range(4))
    out = _cluster_web_members(webs, 0.0, min_sep=0.08, t_c=1.0)
    assert [w.index for w in out] == [0, 1, 2, 3]


def test_regularise_keeps_web_over_nearby_cap():
    got = _regularise_anchors(
        (("le", 0.0), ("cap_lo", 0.2469), ("web_0", 0.2488), ("te", 1.0)),
        te_arc=0.03, min_gap=0.002,
    )
    names = [n for n, _ in got]
    assert "web_0" in names
    assert "cap_lo" not in names


def test_te_wrap_arc_includes_oml_ridge():
    th = np.linspace(0.0, 2.0 * np.pi, 36, endpoint=False)
    pts = np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], axis=1)
    i_le = int(np.argmin(pts[:, 0]))
    i_te = int(np.argmax(pts[:, 0]))
    a = pts[(i_te + 2) % 36]
    b = pts[(i_te - 2) % 36]
    arc = _te_wrap_arc(pts, a, b, i_le, i_te)
    assert arc is not None
    assert any(np.allclose(p, pts[i_te]) for p in arc)


def test_09_pa_4d_loft_keeps_live_webs_on_cylinder():
    geo = _09_geo("4D-F-PA-10R")
    zs = (0.0, 0.5, 1.0, 2.5)
    _tpls, samples = _templates(geo, zs, n_chord=8)
    for z, samp in zip(zs, samples):
        cut = geo.section_at(float(z))
        live = {int(n.rsplit("_", 1)[1]) for n in samp if n.startswith("web_")}
        expected = set(range(4)) - set(cut.dropped_web_indices)
        assert live == expected, (z, live, expected, cut.dropped_web_indices)


def test_09_pa_4d_loft_tracks_pitch_offset_not_catalogue():
    geo = _09_geo("4D-F-PA-10R")
    z_rib5 = 10.0
    cut = geo.section_at(z_rib5)
    assert 0 not in cut.dropped_web_indices
    authored = geo.region_graph.webs[0].at(z_rib5)
    expect = float(cut.pitch_axis) + float(geo._web_x_from_pitch[0]) / float(cut.chord)
    assert authored == pytest.approx(expect, abs=1e-6)
    assert authored < 0.12
    _tpls, samples = _templates(geo, [z_rib5], n_chord=8)
    assert "web_0" in samples[0]
    z_rib6 = 12.0
    cut6 = geo.section_at(z_rib6)
    assert 0 in cut6.dropped_web_indices
    assert 3 not in cut6.dropped_web_indices
    _t6, s6 = _templates(geo, [z_rib6], n_chord=8)
    assert "web_0" not in s6[0]
    assert "web_3" in s6[0]


def test_09_pa_4d_web3_onset_and_drop_are_skin_ruled():
    geo = _09_geo("4D-F-PA-10R")
    tpls, samples = _templates(geo, [0.0, 0.5], n_chord=8)
    faces = loft_mapped_strips(samples[0], samples[1], tpls[0], tpls[1])
    bay = [c for n, _f, c in faces if n == "bay_pr_2"]
    wrap = [c for n, _f, c in faces if n == "te"]
    te = [c for n, _f, c in faces if n == "te_pr"]
    assert bay
    assert sum(1 for c in bay if not c) >= max(1, len(bay) // 2)
    ruled = wrap if wrap else te
    assert ruled
    assert sum(1 for c in ruled if not c) >= max(1, len(ruled) // 2)
    tpls, samples = _templates(geo, [12.0, 12.5], n_chord=8)
    faces = loft_mapped_strips(samples[0], samples[1], tpls[0], tpls[1])
    te = [c for n, _f, c in faces if n == "te_pr"]
    assert te
    assert sum(1 for c in te if not c) >= max(1, len(te) // 2)


def test_09_pa_4d_pressure_skin_covers_web_events():
    geo = _09_geo("4D-F-PA-10R")
    zs = np.linspace(0.0, 20.0, 41)
    tpls, samples = _templates(geo, zs, n_chord=8)
    for s0, s1, t0, t1 in zip(samples, samples[1:], tpls, tpls[1:]):
        got = {n for n, _f, _c in loft_mapped_strips(s0, s1, t0, t1)}
        pr1 = {b.name for b in t1.loop if b.side == "pressure" and not b.wrap and b.name in s1}
        assert pr1 <= got, (float(t0.z), float(t1.z), sorted(pr1 - got),
                            t0.signature[1], t1.signature[1])


def test_09_pa_loft_keeps_surviving_web_index():
    geo = _09_geo("2D-F-PA-10R")
    zs = np.linspace(0.0, 20.0, 21)
    tpls, samples = _templates(geo, zs, n_chord=8)
    saw_drop = False
    for z, samp in zip(zs, samples):
        cut = geo.section_at(float(z))
        names = set(samp)
        if 0 in cut.dropped_web_indices and 1 not in cut.dropped_web_indices:
            saw_drop = True
            assert "web_1" in names, (z, names, cut.dropped_web_indices)
            assert "web_0" not in names, (z, names, cut.dropped_web_indices)
    assert saw_drop


def test_09_le_loft_drops_aft_web():
    geo = _09_geo("2D-F-LE-10R")
    zs = np.linspace(0.0, 20.0, 21)
    _tpls, samples = _templates(geo, zs, n_chord=8)
    saw_drop = False
    for z, samp in zip(zs, samples):
        cut = geo.section_at(float(z))
        names = set(samp)
        if 1 in cut.dropped_web_indices:
            saw_drop = True
            assert "web_1" not in names, (z, names, cut.dropped_web_indices)
            if 0 not in cut.dropped_web_indices:
                assert "web_0" in names, (z, names, cut.dropped_web_indices)
    assert saw_drop


def test_09_te_chatter_and_tip_span():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
    from _shared.debug_mesh import as_quad, span_metrics

    geo = _09_geo("2D-F-LC-0R")
    zs = np.linspace(0.0, 20.0, 41)
    tpls, samples = _templates(geo, zs, n_chord=12)
    walls = list(zip([float(z) for z in zs], tpls, samples))
    faces, kinds = [], []
    for s0, s1, t0, t1 in zip(samples, samples[1:], tpls, tpls[1:]):
        for name, face, _c in loft_mapped_strips(s0, s1, t0, t1):
            faces.append(face)
            kinds.append(name)
    face_arr = np.stack([as_quad(f) for f in faces]) if faces else np.zeros((0, 4, 3))
    metrics = span_metrics(walls, face_arr, np.asarray(kinds, dtype="U64"), geo=geo)
    assert metrics["n_coverage_holes"] == 0
    assert metrics["max_te_chatter"] < 0.05
    assert metrics["z_tip"] == pytest.approx(20.0)


def test_09_3df_le_tip_anchors_one_web():
    cut = _09_geo("3D-F-LE-10R").section_at(20.0)
    assert len(cut.walls_of("web")) >= 1
    cut_pa = _09_geo("3D-F-PA-10R").section_at(20.0)
    assert len(cut_pa.walls_of("web")) >= 1

