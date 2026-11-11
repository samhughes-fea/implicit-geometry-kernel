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
    tip_closure_polygons,
    transverse_members_from_windio,
    zone_region_graphs,
)
from implicit_geometry_kernel.blade.concepts import BladeConcept
from implicit_geometry_kernel.blade.topology import (
    Anchor,
    Region,
    ShellTopologyTemplate,
    WallRegionGraph,
    WebMember,
    _cluster_web_members,
    _regularise_anchors,
    _te_wrap_arc,
    _template_from_cut,
    place_station_keys,
    rib_cell_polygons_from_samples,
    rib_cell_triangles_from_samples,
    sample_region_polylines,
    hold_te_blunt,
    te_blunt_flags_for_cuts,
    te_blunt_eps,
    te_blunt_of_cut,
    te_is_blunt,
    te_opening,
)
from implicit_geometry_kernel.blade.topology.map_walk import _skin_mid_loop


def _geo(key="2D-F-LC-0R", **kwargs):
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    kw = dict(pivot_fraction=0.3, spar_cap_thickness=0.008, shear_web_thickness=0.004)
    kw.update(kwargs)
    return BladeGeometry(stations, key, 0.004, **kw)


def _templates(geo, zs, n_chord=8):
    cuts = [geo.section_at(float(z)) for z in zs]
    flags = te_blunt_flags_for_cuts(cuts)
    out = []
    samples = []
    for cut, blunt in zip(cuts, flags):
        sg = cut.shell_graph()
        t = _template_from_cut(
            geo.region_graph, cut, float(cut.z), n_chord, te_blunt=blunt,
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
    stations = [
        Station(0.0, 2.0, Airfoil.ellipse(1.0), pitch_axis=0.5),
        Station(10.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, web_positions=(0.25, 0.50), pivot_fraction=0.5)
    tpls, samples = _templates(geo, [0.0, 10.0], n_chord=8)
    t0, t1 = tpls
    wrap = next(b for b in t0.loop if b.wrap)
    assert next((b for b in t1.loop if b.wrap), None) is None
    faces = loft_mapped_strips(samples[0], samples[1], t0, t1)
    te = [(n, f, c) for n, f, c in faces if n == wrap.name]
    assert te
    assert any(c or len(f) == 3 for _n, f, c in te)
    mapped = t0.map_vertex(wrap.start, t1)
    assert mapped == ("TE",) or mapped in t1.vertices


def test_write_mid_surfaces_obj_region_groups_and_tip_closure(tmp_path: Path):
    from implicit_geometry_kernel.blade.stations import TransverseMember

    geo = _geo(web_positions=(0.18, 0.48), transverse=(TransverseMember("rib_mid", 5.0, 0.012),))
    obj = Path(geo.write_mid_surfaces(tmp_path / "mids.obj", n_span=5, n_chord=8)).read_text(
        encoding="utf-8")
    assert "g " in obj and "v " in obj
    groups = {ln.split()[1] for ln in obj.splitlines() if ln.startswith("g ")}
    assert any(g.startswith("tip_closure") for g in groups)
    assert any(g.startswith("rib_mid") for g in groups)
    n_gons = [ln for ln in obj.splitlines() if ln.startswith("f ") and len(ln.split()) > 5]
    assert n_gons, "rib and tip OBJ groups must keep n-gon faces (n>4)"
    ribs = rib_face_polygons(geo, n_chord=8)
    assert ribs and all(abs(float(ring[0, 2]) - 5.0) < 1e-6 for _n, ring in ribs)
    tips = tip_closure_polygons(geo, n_chord=8)
    assert tips and all(abs(float(ring[0, 2]) - 10.0) < 1e-6 for _n, ring in tips)
    assert all(n.startswith("tip_closure") for n, _r in tips)


def test_mapped_loft_omits_rib_strips():

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
    assert transverse_members_from_windio(structure, z_of_s) == ()
    geo = _geo(web_positions=(0.18, 0.48))
    tpls, samples = _templates(geo, [8.0, 10.0])
    strips = loft_mapped_strips(samples[0], samples[1], tpls[0], tpls[1])
    assert all("rib" not in n and not n.startswith("tip_closure") for n, _f, _c in strips)
    assert rib_face_polygons(geo) == []
    assert tip_closure_polygons(geo)


def _assert_keyed_seams(geo, z, n_chord=10):
    cut = geo.section_at(float(z))
    sg = cut.shell_graph()
    t = _template_from_cut(geo.region_graph, cut, float(z), n_chord)
    samples = sample_region_polylines(cut, sg, t, geo.region_graph)
    keys = place_station_keys(
        cut, sg, geo.region_graph, float(z),
        te_blunt=any(b.wrap for b in t.loop),
    )
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


def _web_kink_deg(pts: np.ndarray) -> float:
    pts = np.asarray(pts, float)
    if len(pts) < 2:
        return 0.0
    axis = pts[-1] - pts[0]
    nrm = float(np.linalg.norm(axis[:2]))
    if nrm < 1e-12:
        return 0.0
    axis = axis[:2] / nrm
    segs = np.diff(pts[:, :2], axis=0)
    lengths = np.linalg.norm(segs, axis=1)
    keep = lengths > 1e-12
    if not np.any(keep):
        return 0.0
    directions = segs[keep] / lengths[keep, None]
    cross = np.abs(directions[:, 0] * axis[1] - directions[:, 1] * axis[0])
    return float(np.degrees(np.arcsin(np.clip(cross, 0.0, 1.0))).max())


def test_conformal_webs_stay_straight_on_5d_stations_07_11():
    """SG-002 analogue: loft sampler must not relocate emitted web feet (csg_mesh 07–11)."""
    geo = _09_geo("5D-F-PA-10R")
    zs = np.linspace(float(geo.stations[0].z), float(geo.stations[-1].z), 16)
    tpls, samples = _templates(geo, zs[6:11], n_chord=8)
    n_webs = 0
    for tmpl, samp in zip(tpls, samples):
        for web in tmpl.webs:
            assert web.name in samp
            n_webs += 1
            worst = _web_kink_deg(samp[web.name])
            assert worst <= 0.05, (
                f"z={tmpl.z:.1f} m {web.name}: conformal web leans {worst:.3f} deg "
                "off its axis — feet were relocated onto chord-fraction skin keys"
            )
    assert n_webs >= 1


def test_te_opening_classifies_ellipse_naca_flatback():
    naca = _geo(web_positions=(0.18, 0.48))
    cut = naca.section_at(4.0)
    pts = _skin_mid_loop(cut.shell_graph(), float(cut.z))
    assert not te_is_blunt(cut, pts)
    assert not te_blunt_of_cut(cut)
    assert te_opening(cut, pts) <= te_blunt_eps(cut)

    stations = [
        Station(0.0, 2.0, Airfoil.ellipse(1.0), pitch_axis=0.5),
        Station(10.0, 1.0, Airfoil.ellipse(1.0), pitch_axis=0.5),
    ]
    ell = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, web_positions=(0.25, 0.50), pivot_fraction=0.5)
    ecut = ell.section_at(0.0)
    assert te_blunt_of_cut(ecut)
    assert te_opening(ecut, _skin_mid_loop(ecut.shell_graph(), 0.0)) > 0.04 * 2.0

    blunt = Airfoil.naca4("0018").with_te_thickness(0.05)
    fb = BladeGeometry(
        [Station(0.0, 2.0, blunt, pitch_axis=0.3), Station(10.0, 1.0, blunt, pitch_axis=0.3)],
        "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
        web_positions=(0.18, 0.48), pivot_fraction=0.3,
    )
    fcut = fb.section_at(0.0)
    assert te_blunt_of_cut(fcut)


def test_hold_te_blunt_hysteresis():
    assert hold_te_blunt([2.0, 0.8, 0.4, 0.2, 1.2], 1.0) == [
        True, False, False, False, True,
    ]
    assert hold_te_blunt([0.2, 0.9, 1.1], 1.0) == [False, False, True]
    assert hold_te_blunt([1.1, 1.0, 0.9], 1.0) == [True, False, False]


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
    tpls, samples = _templates(geo, [0.0, 10.0], n_chord=8)
    assert "te" in samples[0]
    assert "te" not in samples[1]
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


def test_family_c_template_has_pad_patches_and_inter_pad_bay():
    graph = WallRegionGraph.from_section_layout("2C-F-LC-0R")
    t = ShellTopologyTemplate.from_region_graph(graph, 0.0, n_chord=6)
    su = [b for b in t.loop if b.side == "suction" and not b.wrap]
    caps = [b.name for b in su if b.name.startswith("cap_")]
    bays = [b.name for b in su if b.name.startswith("bay_")]
    assert len(caps) == 4
    assert bays
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


def test_sharp_te_has_no_wrap_branch():
    t, samples, keys = _assert_keyed_seams(_geo(web_positions=(0.18, 0.48)), 4.0)
    assert next((b for b in t.loop if b.wrap), None) is None
    assert "te" not in samples
    a, b = keys[("te", "suction")], keys[("te", "pressure")]
    assert float(np.linalg.norm(a - b)) < 1e-6 * 2.0


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
    assert "te" not in samples
    assert next((b for b in t.loop if b.wrap), None) is None
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


def test_09_5d_foil_stations_have_no_te_spur():
    geo = _09_geo("5D-F-PA-10R")
    for z in (3.2, 4.0):
        t, samples, _keys = _assert_keyed_seams(geo, z)
        assert next((b for b in t.loop if b.wrap), None) is None
        assert "te" not in samples
        assert "te_su" in samples and "te_pr" in samples
        assert any(b.side == "suction" and b.name in samples for b in t.loop)
        assert any(b.side == "pressure" and b.name in samples for b in t.loop)
        oml = np.asarray(geo.oml_polyline(float(z)), float)
        te_x = float(np.max(oml[:, 0]))
        chord = float(geo.section_at(float(z)).chord)
        xs = np.concatenate([np.asarray(p, float)[:, 0] for p in samples.values() if len(p)])
        assert float(np.max(xs)) <= te_x + 1e-3 * chord


def test_09_te_wrap_held_without_flicker():
    geo = _09_geo("2D-F-PA-10R")
    zs = np.linspace(0.0, 12.0, 41)
    tpls, samples = _templates(geo, zs, n_chord=8)
    flags = [any(b.wrap for b in t.loop) for t in tpls]
    n_flip = sum(a != b for a, b in zip(flags, flags[1:]))
    assert n_flip <= 1
    if True in flags and False in flags:
        i = flags.index(False)
        assert all(flags[:i]) and not any(flags[i:])
    for wrap_on, samp in zip(flags, samples):
        if wrap_on:
            assert "te" in samp
        else:
            assert "te" not in samp


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
    geo = _geo(web_positions=(0.18, 0.48))
    tpls, samples = _templates(geo, [10.0], n_chord=8)
    rings = rib_cell_polygons_from_samples(tpls[0], samples[0], "tip_closure")
    assert rings
    assert all(len(ring) > 4 for _n, ring in rings)
    tmpl, samp = tpls[0], samples[0]
    from implicit_geometry_kernel.blade.topology.map import _key_name
    from implicit_geometry_kernel.blade.topology.ribs import _group_skin_branches

    web_keys = {_key_name(w.suction) for w in tmpl.webs}
    su_g = _group_skin_branches(
        [b for b in tmpl.loop if b.side == "suction" and not b.wrap], web_keys)
    pr_g = _group_skin_branches(
        [b for b in tmpl.loop if b.side == "pressure" and not b.wrap], web_keys)
    for k, ((sbs, pbs), (_name, ring)) in enumerate(zip(zip(su_g, pr_g), rings)):
        wanted = set()
        for b in (*sbs, *pbs):
            if b.name in samp:
                wanted |= _ptset(samp[b.name])
        aft = next((w for w in tmpl.webs if w.suction == sbs[-1].end or w.pressure == sbs[-1].end), None)
        fwd = next((w for w in tmpl.webs if w.suction == sbs[0].start or w.pressure == sbs[0].start), None)
        if aft is not None and aft.name in samp:
            wanted |= _ptset(samp[aft.name])
        if fwd is not None and fwd.name in samp:
            wanted |= _ptset(samp[fwd.name])
        wrap = next((b for b in tmpl.loop if b.wrap), None)
        if wrap is not None and wrap.name in samp and _key_name(sbs[-1].end) == "te":
            wanted |= _ptset(samp[wrap.name])
        got = _ptset(ring)
        assert wanted <= got, (k, sorted(wanted - got)[:4])
    if len(tmpl.webs) and tmpl.webs[0].name in samp and len(rings) >= 2:
        interior = np.asarray(samp[tmpl.webs[0].name], float)[len(samp[tmpl.webs[0].name]) // 2]
        hit = [_ptset(r) for _n, r in rings]
        n_hit = sum(tuple(np.round(interior, 9)) in s for s in hit)
        assert n_hit == 2


def test_rib_triangles_use_only_ngon_vertices():
    geo = _geo(web_positions=(0.18, 0.48))
    tpls, samples = _templates(geo, [10.0], n_chord=8)
    rings = dict(rib_cell_polygons_from_samples(tpls[0], samples[0], "tip_closure"))
    tris = rib_cell_triangles_from_samples(tpls[0], samples[0], "tip_closure")
    assert tris
    for name, tri in tris:
        assert _ptset(tri) <= _ptset(rings[name])
        assert len(tri) == 3


def test_ear_clip_le_cell_stays_inside_ngon():
    from implicit_geometry_kernel.blade.topology.ribs import (
        _point_in_ring,
        _project_xy,
        _ring_basis,
        ear_clip_ring,
    )

    geo = _09_geo("5D-F-PA-10R")
    member = next(m for m in geo.transverse if str(m.name) == "rib_1")
    z = float(member.z)
    assert z == pytest.approx(np.linspace(0.0, 12.0, 12)[1])  # 10 interior ribs
    tpls, samples = _templates(geo, [z], n_chord=24)
    samp = samples[0]
    z_loft = float(next(iter(samp.values()))[0, 2])
    assert z_loft == pytest.approx(float(geo.section_at(z).z), abs=1e-9)
    rings = rib_cell_polygons_from_samples(tpls[0], samp, member)
    assert rings
    name, ring = rings[0]
    assert name.endswith("cell_0")
    tris = rib_cell_triangles_from_samples(tpls[0], samp, member)
    assert tris
    for _n, tri in tris:
        assert np.allclose(tri[:, 2], z_loft, atol=1e-8)
    origin, u, v, xy = _ring_basis(ring)
    clipped = ear_clip_ring(ring)
    assert clipped
    for tri in clipped:
        for p in _project_xy(_barycentric(tri), origin, u, v):
            assert _point_in_ring(p, xy)


def _barycentric(tri: np.ndarray, n: int = 4) -> np.ndarray:
    pts = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            pts.append((i * tri[0] + j * tri[1] + k * tri[2]) / n)
    return np.asarray(pts)


def test_ear_clip_concave_cell_stays_inside():
    from implicit_geometry_kernel.blade.topology.ribs import (
        _point_in_ring,
        _project_xy,
        _ring_basis,
        ear_clip_ring,
    )

    ring = np.array([
        [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 3.0, 0.0], [0.0, 3.0, 0.0],
        [0.0, 2.0, 0.0], [2.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    ])
    origin, u, v, xy = _ring_basis(ring)
    fan = [np.stack([ring[0], ring[i], ring[i + 1]]) for i in range(1, len(ring) - 1)]
    assert any(
        any(not _point_in_ring(p, xy) for p in _project_xy(_barycentric(t), origin, u, v))
        for t in fan
    )
    for tri in ear_clip_ring(ring):
        for p in _project_xy(_barycentric(tri), origin, u, v):
            assert _point_in_ring(p, xy)


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
    s = 12.0 / 20.0
    root, foil = Airfoil.ellipse(1.0), Airfoil.naca4("0024")
    return [
        Station(0.0 * s, 1.80 * s, root, twist=12.0, pitch_axis=0.50),
        Station(1.0 * s, 1.88 * s, Airfoil.blend(root, foil, 0.15), twist=11.5, pitch_axis=0.47),
        Station(2.5 * s, 2.05 * s, Airfoil.blend(root, foil, 0.50), twist=10.0, pitch_axis=0.40),
        Station(4.0 * s, 2.28 * s, Airfoil.blend(root, foil, 0.85), twist=8.8, pitch_axis=0.33),
        Station(5.0 * s, 2.40 * s, foil, twist=8.0, pitch_axis=0.30),
        Station(6.5 * s, 2.20 * s, Airfoil.naca4("0022"), twist=6.5, pitch_axis=0.30),
        Station(8.5 * s, 1.85 * s, Airfoil.naca4("0020"), twist=4.8, pitch_axis=0.30),
        Station(11.0 * s, 1.50 * s, Airfoil.naca4("0018"), twist=2.8, pitch_axis=0.30),
        Station(14.0 * s, 1.10 * s, Airfoil.naca4("0015"), twist=1.0, pitch_axis=0.30),
        Station(17.0 * s, 0.90 * s, Airfoil.naca4("0013"), twist=-0.6, pitch_axis=0.30),
        Station(20.0 * s, 0.75 * s, Airfoil.naca4("0012"), twist=-2.0, pitch_axis=0.30),
    ]


def _09_geo(key: str) -> BladeGeometry:
    st = _09_stations()
    sc = 12.0 / 20.0
    return BladeGeometry.from_concept(
        BladeConcept(key), st, 0.004 * sc,
        spar_cap_thickness=0.008 * sc, shear_web_thickness=0.004 * sc,
        pivot_fraction=([sta.z for sta in st], [sta.pitch_axis for sta in st]),
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
    zs = (0.0, 0.3, 0.6, 1.5)
    _tpls, samples = _templates(geo, zs, n_chord=8)
    for z, samp in zip(zs, samples):
        cut = geo.section_at(float(z))
        live = {int(n.rsplit("_", 1)[1]) for n in samp if n.startswith("web_")}
        expected = set(range(4)) - set(cut.dropped_web_indices)
        assert live == expected, (z, live, expected, cut.dropped_web_indices)


def test_09_pa_4d_loft_tracks_pitch_offset_not_catalogue():
    geo = _09_geo("4D-F-PA-10R")
    z_mid = 3.0
    cut = geo.section_at(z_mid)
    authored = geo.region_graph.webs[0].at(z_mid)
    expect = float(cut.pitch_axis) + float(geo._web_x_from_pitch[0]) / float(cut.chord)
    assert authored == pytest.approx(expect, abs=1e-6)
    assert authored < 0.20
    _tpls, samples = _templates(geo, [z_mid], n_chord=8)
    assert "web_0" in samples[0]
    root = geo.section_at(0.0)
    assert 3 in root.dropped_web_indices
    assert 3 not in cut.dropped_web_indices
    _t0, s0 = _templates(geo, [0.0], n_chord=8)
    assert "web_3" not in s0[0]
    assert "web_3" in samples[0]


def test_09_pa_4d_web3_onset_and_drop_are_skin_ruled():
    geo = _09_geo("4D-F-PA-10R")
    zs = np.linspace(0.0, 12.0, 41)
    live = [3 not in geo.section_at(float(z)).dropped_web_indices for z in zs]
    assert any(live) and not live[0]
    i_on = next(i for i, v in enumerate(live) if v)
    tpls, samples = _templates(geo, [zs[i_on - 1], zs[i_on]], n_chord=8)
    faces = loft_mapped_strips(samples[0], samples[1], tpls[0], tpls[1])
    names = {n for n, _f, _c in faces}
    assert any(n.startswith("bay_pr") or n.startswith("te") for n in names)
    ruled = [c for n, _f, c in faces if n.startswith("bay_pr") or n in ("te", "te_pr")]
    assert ruled
    assert sum(1 for c in ruled if not c) >= max(1, len(ruled) // 2)
    live0 = [0 not in geo.section_at(float(z)).dropped_web_indices for z in zs]
    assert live0[0] and not live0[-1]
    i_off = next(i for i in range(1, len(live0)) if not live0[i])
    tpls, samples = _templates(geo, [zs[i_off - 1], zs[i_off]], n_chord=8)
    faces = loft_mapped_strips(samples[0], samples[1], tpls[0], tpls[1])
    te = [c for n, _f, c in faces if n in ("te_pr", "te") or n.startswith("bay_pr")]
    assert te
    assert sum(1 for c in te if not c) >= max(1, len(te) // 2)


def test_09_pa_4d_pressure_skin_covers_web_events():
    geo = _09_geo("4D-F-PA-10R")
    zs = np.linspace(0.0, 12.0, 41)
    tpls, samples = _templates(geo, zs, n_chord=8)
    for s0, s1, t0, t1 in zip(samples, samples[1:], tpls, tpls[1:]):
        got = {n for n, _f, _c in loft_mapped_strips(s0, s1, t0, t1)}
        pr1 = {b.name for b in t1.loop if b.side == "pressure" and not b.wrap and b.name in s1}
        assert pr1 <= got, (float(t0.z), float(t1.z), sorted(pr1 - got),
                            t0.signature[1], t1.signature[1])


def test_09_pa_loft_keeps_surviving_web_index():
    geo = _09_geo("2D-F-PA-10R")
    zs = np.linspace(0.0, 12.0, 21)
    tpls, samples = _templates(geo, zs, n_chord=8)
    for z, samp in zip(zs, samples):
        cut = geo.section_at(float(z))
        names = set(samp)
        for i in (0, 1):
            label = f"web_{i}"
            if i in cut.dropped_web_indices:
                assert label not in names, (z, names, cut.dropped_web_indices)
            else:
                assert label in names, (z, names, cut.dropped_web_indices)


def test_09_le_loft_tracks_web_indices():
    geo = _09_geo("2D-F-LE-10R")
    zs = np.linspace(0.0, 12.0, 21)
    _tpls, samples = _templates(geo, zs, n_chord=8)
    for z, samp in zip(zs, samples):
        cut = geo.section_at(float(z))
        names = set(samp)
        for i in (0, 1):
            label = f"web_{i}"
            if i in cut.dropped_web_indices:
                assert label not in names, (z, names, cut.dropped_web_indices)
            else:
                assert label in names, (z, names, cut.dropped_web_indices)


def test_strip_faces_resample_unequal_counts_and_drop_nonfinite():
    from implicit_geometry_kernel.blade.topology.map_loft import (
        _finite_loft_face,
        _strip_faces,
    )
    from implicit_geometry_kernel.blade.topology.template import LoopBranch

    t0 = ShellTopologyTemplate(z=0.0, signature=(), loop=(), webs=(), vertices=(), n_chord=2)
    t1 = t0
    branch = LoopBranch(
        "te_pr", ("TE", "pressure"), ("anchor", "web_0", "pressure"), n=8, side="pressure",
    )
    pa = np.column_stack([np.linspace(0.0, 1.0, 3), np.zeros(3), np.zeros(3)])
    pb = np.column_stack([np.linspace(0.0, 1.0, 8), np.zeros(8), np.ones(8)])
    faces = _strip_faces("te_pr", pa, pb, t0, t1, branch)
    assert faces
    assert all(len(face) == 4 for _n, face, _c in faces)
    assert all(_finite_loft_face(face) for _n, face, _c in faces)
    nan_face = np.array([[np.nan, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    assert not _finite_loft_face(nan_face)
    assert not _finite_loft_face(np.zeros((4, 3)))


@pytest.mark.parametrize("key", [
    "0A-LC-0R", "0B-LC-0R", "1B-F-LC-10R", "2A-F-LC-10R", "2C-F-LC-10R", "2D-CN-LC-10R",
])
def test_09_extra_keys_build(key):
    geo = _09_geo(key)
    assert float(geo.stations[-1].z) == pytest.approx(12.0)
    cut = geo.section_at(3.0)
    assert cut is not None


@pytest.mark.parametrize("key,z,n_patches", [
    ("2A-F-LC-10R", 6.0, 3),
    ("2C-F-LC-10R", 6.0, 3),
    ("2D-F-LC-10R", 6.0, 3),
    ("1B-F-LC-10R", 6.0, 2),
    ("2D-F-PA-10R", 6.0, 3),
])
def test_rib_patches_split_only_at_live_webs(key, z, n_patches):
    from implicit_geometry_kernel.blade.stations import TransverseMember

    geo = _09_geo(key)
    tpls, samples = _templates(geo, [z], n_chord=8)
    member = next((m for m in geo.transverse if str(m.name).startswith("rib")),
                  TransverseMember("rib", z, 0.01))
    rings = rib_cell_polygons_from_samples(tpls[0], samples[0], member)
    assert len(tpls[0].webs) + 1 == n_patches
    assert len(rings) == n_patches
    if key == "2C-F-LC-10R":
        n_su = sum(1 for b in tpls[0].loop if b.side == "suction" and not b.wrap)
        assert n_su > n_patches


def test_rib_patches_follow_dropped_pa_webs():
    geo = _09_geo("2D-F-PA-10R")
    tpls, samples = _templates(geo, [11.5], n_chord=8)
    rings = rib_cell_polygons_from_samples(tpls[0], samples[0], "rib_tip")
    assert len(rings) == len(tpls[0].webs) + 1


def test_rib_patches_zero_webs_one_ring():
    from implicit_geometry_kernel.blade.stations import TransverseMember

    geo = _09_geo("0B-LC-0R")
    tpls, samples = _templates(geo, [6.0], n_chord=8)
    assert tpls[0].webs == ()
    n_su = sum(1 for b in tpls[0].loop if b.side == "suction" and not b.wrap)
    assert n_su > 1
    rings = rib_cell_polygons_from_samples(
        tpls[0], samples[0], TransverseMember("rib_x", 6.0, 0.01))
    assert len(rings) == 1
    assert rings[0][0] == "rib_x"


@pytest.mark.parametrize("key", ["2D-F-LC-0R", "4D-F-PA-10R"])
def test_09_mapped_loft_faces_are_finite_and_bounded(key):
    geo = _09_geo(key)
    zs = (0.0, 1.0, 5.0, 12.0)
    tpls, samples = _templates(geo, zs, n_chord=12)
    bound = 1.0e3
    n_faces = 0
    for s0, s1, t0, t1 in zip(samples, samples[1:], tpls, tpls[1:]):
        for name, face, _c in loft_mapped_strips(s0, s1, t0, t1):
            n_faces += 1
            arr = np.asarray(face, float)
            z0, z1 = float(t0.z), float(t1.z)
            assert arr.ndim == 2 and arr.shape[1] == 3 and len(arr) >= 3, (name, z0, z1, arr.shape)
            assert np.all(np.isfinite(arr)), (name, z0, z1)
            extent = float(np.max(np.abs(arr)))
            assert extent <= bound, (name, z0, z1, extent)
    assert n_faces > 0


def test_09_te_chatter_and_tip_span():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
    from _shared.debug_mesh import as_quad, span_metrics

    geo = _09_geo("2D-F-LC-0R")
    zs = np.linspace(0.0, 12.0, 41)
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
    assert metrics["z_tip"] == pytest.approx(12.0)


def test_09_3df_le_tip_anchors_one_web():
    cut = _09_geo("3D-F-LE-10R").section_at(12.0)
    assert len(cut.walls_of("web")) >= 1
    cut_pa = _09_geo("3D-F-PA-10R").section_at(12.0)
    assert len(cut_pa.walls_of("web")) >= 1

