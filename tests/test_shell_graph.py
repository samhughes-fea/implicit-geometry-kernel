"""Two-view rule, loftable shell graph, merged wall_meta, stack-aware wall_mid."""

import numpy as np
import pytest

from implicit_geometry_kernel import (
    Airfoil,
    BladeGeometry,
    Station,
    Stringer,
    SYSTEM_TYPE_LAYOUT_KEYS,
    TransverseMember,
    partition_edge_curve,
    shared_loft_names,
    station_network,
    topology_diff,
    wall_meta_from_cut,
)
from implicit_geometry_kernel.blade.shell_graph import ShellGraph, ShellRun
from implicit_geometry_kernel.contract import InfeasibleSectionGeometry
from implicit_geometry_kernel.sdf import FeatureTable


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                         shear_web_thickness=0.002, web_positions=(0.18, 0.48), pivot_fraction=0.3)


def assert_closed_conformal_loop(cut, *, cap_on_loop: bool = True):
    g = cut.shell_graph()
    names = [r.name for r in g.loop]
    assert names, "empty skin loop"
    assert len(names) == len(set(names)), names
    web_names = [r.name for r in g.webs]
    assert len(web_names) == len(set(web_names)), web_names
    all_names = [r.name for r in g.runs()]
    assert len(all_names) == len(set(all_names)), all_names
    chain = g.closed_loop_chain()
    assert chain is not None, names
    if not g.loop:
        return
    scale = max(1.0, max(float(np.max(np.abs(r.points[:, :2]))) for r in g.loop))
    tol = 1e-6 * scale
    loop_xy = np.vstack([np.asarray(r.points[:, :2], float) for r in g.loop])
    for web in g.webs:
        for i in (0, -1):
            d = float(np.min(np.linalg.norm(loop_xy - web.points[i, :2], axis=1)))
            assert d <= tol, (web.name, i, d)
    skins = list(cut.walls_of("skin"))
    if skins:
        skin_xy = np.vstack([w.xy for w in skins])
        for run in g.loop:
            if run.kind != "cap":
                continue
            d = np.linalg.norm(
                run.points[:, None, :2] - skin_xy[None, :, :], axis=2
            ).min(axis=1)
            assert float(np.max(d)) <= max(tol, 1e-4), (run.name, float(np.max(d)))
    if cap_on_loop:
        assert any(r.kind == "cap" for r in g.loop)


def _assert_unique_closed_loop(cut, *, cap_on_loop: bool = True):
    assert_closed_conformal_loop(cut, cap_on_loop=cap_on_loop)


def test_walls_keep_caps_mid_surface_does_not():
    cut = _geo().section_at(4.0)
    assert any(w.kind == "cap" for w in cut.walls)
    ms = cut.mid_surface()
    assert {s.kind for s in ms.segments} <= {"skin", "web"}
    assert "cap" not in {s.kind for s in ms.segments}


def test_two_web_loop_run_names_unique_and_closed():
    _assert_unique_closed_loop(_geo().section_at(4.0))


def test_cap_windows_span_published_chord_interval():
    """Published (lo, hi) is the cap band, not a pitch-centred pad of half the width."""
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
        shear_web_thickness=0.002, web_positions=(0.18, 0.48),
        cap_windows=((0.20, 0.50),), pivot_fraction=0.3)
    cut = geo.section_at(0.0)
    _assert_unique_closed_loop(cut)
    caps = [r for r in cut.shell_graph().loop if r.kind == "cap"]
    assert caps
    window = 0.30 * cut.chord
    for run in caps:
        span = float(np.ptp(np.asarray(run.points[:, :2], float)[:, 0]))
        assert span >= 0.5 * window, (run.name, span, window)
    solids = list(cut.walls_of("cap"))
    assert solids
    for wall in solids:
        span = float(np.ptp(wall.xy[:, 0]))
        assert span >= 0.5 * window, (wall.label, span, window)


def test_one_web_loop_run_names_unique_and_closed():
    """Cap inside :leading_edge must not emit two loop runs with one name."""
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "1C-F-LC-0R", 0.003, spar_cap_thickness=0.008,
        shear_web_thickness=0.002, web_positions=(0.30,), pivot_fraction=0.3)
    cut = geo.section_at(4.0)
    assert sum(1 for w in cut.walls if w.kind == "web") == 1
    _assert_unique_closed_loop(cut)


def test_one_web_chord_window_keeps_unindexed_cap_name():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "1C-F-LC-0R", 0.003, spar_cap_thickness=0.008,
        shear_web_thickness=0.002, web_positions=(0.30,),
        cap_windows=((0.20, 0.50),), pivot_fraction=0.3)
    cut = geo.section_at(4.0)
    assert_closed_conformal_loop(cut)
    names = [r.name for r in cut.shell_graph().loop]
    cap_names = [n for n in names if n.startswith("spar_cap_")]
    assert cap_names
    assert not any(n.startswith("spar_cap_suction_surface_0") or n.startswith("spar_cap_pressure_surface_0")
                   for n in cap_names)


def test_shell_graph_loop_partitions_and_combined_thickness():
    cut = _geo().section_at(4.0)
    g = cut.shell_graph()
    assert {r.surface for r in g.loop} == {"skin"}
    assert {r.kind for r in g.loop} <= {"skin", "cap"}
    assert any(r.kind == "cap" for r in g.loop)
    assert all(r.kind == "web" for r in g.webs)
    assert all(r.surface == "web" for r in g.webs)
    cap = next(r for r in g.loop if r.kind == "cap")
    host = next(lb for lb in cap.layers if lb.startswith("skin"))
    cap_label = next(lb for lb in cap.layers if not lb.startswith("skin"))
    t_cap = cut.wall(cap_label).thickness
    t_skin = cut.wall(host).thickness
    assert cap.thickness == t_skin
    assert cap.requires_eccentricity
    assert cap.eccentricity_m == pytest.approx(0.5 * t_cap)
    assert cap.layers == (cap_label, host)
    regions = g.offset_regions()
    assert {r.name for r in regions} == {r.name for r in g.loop if r.kind == "cap"}
    assert all(r.requires_eccentricity for r in regions)
    skin = next(r for r in g.loop if r.kind == "skin")
    assert skin.layers == (skin.name,)
    assert skin.thickness == cut.wall(skin.name).thickness
    # no extra closed cap loop: cap colour only on the skin loop
    assert not any(r.surface == "cap" for r in g.runs())
    assert all(r.points.shape[1] == 3 for r in g.runs())
    assert all(np.allclose(r.points[:, 2], cut.z) for r in g.runs())


def test_wall_meta_from_shell_graph_not_nested_cap():
    cut = _geo().section_at(4.0)
    field, table = wall_meta_from_cut(cut)
    g = cut.shell_graph()
    web = g.webs[0]
    sample = field.evaluate_plane(web.points[len(web.points) // 2, :2].reshape(1, 2), cut.z)
    assert sample.phi[0] < 0.0
    assert table[int(sample.meta[0])].kind == "web"
    assert table[int(sample.meta[0])].name == web.name
    cap = next(r for r in g.loop if r.kind == "cap")
    sample = field.evaluate_plane(cap.points[len(cap.points) // 2, :2].reshape(1, 2), cut.z)
    spec = table[int(sample.meta[0])]
    assert spec.kind == "cap"
    assert spec.name == cap.name
    net = station_network(cut, FeatureTable())
    assert {w.kind for w in net.walls} <= {"skin", "web", "cap"}
    assert any(w.kind == "cap" for w in net.walls)
    assert any(w.kind == "web" for w in net.walls)
    names = [w.label for w in net.walls]
    assert names.count(cap.name) == 1


def test_feature_ids_reuse_shared_run_names():
    geo = _geo()
    table = FeatureTable()
    n0 = station_network(geo.section_at(2.0), table)
    n1 = station_network(geo.section_at(8.0), table)
    shared = {w.label for w in n0.walls} & {w.label for w in n1.walls}
    assert shared
    for name in shared:
        assert table.by_name(name).id == next(w.feature_id for w in n0.walls if w.label == name)
        assert table.by_name(name).id == next(w.feature_id for w in n1.walls if w.label == name)


def test_wall_mid_near_stacked_cap_and_bare_skin():
    geo = _geo()
    cut = geo.section_at(4.0)
    g = cut.shell_graph()
    mid = geo.wall_mid_sdf()
    cap = next(r for r in g.loop if r.kind == "cap")
    p_cap = cap.points[len(cap.points) // 2]
    # Pass 1 bound is run interiors only (not cap/skin joint vertices).
    assert abs(mid.evaluate(p_cap.reshape(1, 3)).phi[0]) < 5e-4
    skin = next(r for r in g.loop if r.kind == "skin")
    p_skin = skin.points[len(skin.points) // 2]
    assert abs(mid.evaluate(p_skin.reshape(1, 3)).phi[0]) < 5e-4
    cavity = geo.cavity_sdf().evaluate(p_cap.reshape(1, 3)).phi[0]
    oml = geo.oml_sdf().evaluate(p_cap.reshape(1, 3)).phi[0]
    assert oml < 0.0
    # Mesh sits on the skin mid; cavity uses stack T, so cavity + oml is t_cap.
    assert cavity + oml == pytest.approx(2.0 * cap.eccentricity_m, abs=1e-3)
    hit = geo.snap_wall_mid((p_cap + np.array([0.0, 0.001, 0.0])).reshape(1, 3))
    assert hit.converged[0]
    np.testing.assert_allclose(hit.points[0, 2], cut.z, atol=1e-12)


def test_partition_edges_include_cap_feature_lines():
    cut = _geo().section_at(4.0)
    g = cut.shell_graph()
    edges = g.partition_edges()
    assert edges
    assert all(e.z == cut.z for e in edges)
    assert all(e.point.shape == (3,) for e in edges)
    caps = {r.name for r in g.loop if r.kind == "cap"}
    assert caps
    named = {(e.name_a, e.name_b) for e in edges}
    assert any(caps & {a, b} for a, b in named)
    cap_name = next(iter(caps))
    hit = g.edge(cap_name, next(e.name_b if e.name_a == cap_name else e.name_a
                                for e in edges if cap_name in (e.name_a, e.name_b)))
    assert hit is not None
    assert np.allclose(hit.point[2], cut.z)


def test_partition_edge_curve_reuses_names_across_span():
    geo = _geo()
    g = geo.section_at(4.0).shell_graph()
    cap = next(e for e in g.partition_edges() if "spar_cap" in e.name_a or "spar_cap" in e.name_b)
    curve = partition_edge_curve(geo, cap.name_a, cap.name_b, [2.0, 4.0, 8.0])
    assert curve.shape == (3, 3)
    np.testing.assert_allclose(curve[:, 2], [2.0, 4.0, 8.0])
    missing = partition_edge_curve(geo, cap.name_a, "not_a_run", [4.0])
    assert missing.shape == (0, 3)


def _run(name, kind, surface="skin"):
    pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    return ShellRun(name, surface, kind, (name,), pts, 0.01)


def test_topology_diff_cap_runout_dropped_web_bay_split():
    a = ShellGraph(0.0, (_run("skin_bay", "skin"), _run("spar_cap", "cap")),
                   (_run("shear_web_0", "web", "web"),))
    b = ShellGraph(1.0, (_run("skin_bay", "skin"), _run("skin_bay_new", "skin")), ())
    diff = topology_diff(a, b)
    assert diff.shared == ("skin_bay",)
    kinds = {e.kind for e in diff.events}
    assert {"cap_runout", "dropped_web", "bay_split"} <= kinds
    assert "spar_cap" not in shared_loft_names([r.name for r in a.runs()], [r.name for r in b.runs()])


def test_topology_diff_dropped_web_on_tapering_blade():
    stations = [Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
                Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3)]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
                        web_offsets=(0.4, 1.2), pivot_fraction=0.3)
    root, tip = geo.section_at(0.0).shell_graph(), geo.section_at(10.0).shell_graph()
    diff = topology_diff(root, tip)
    assert any(e.kind == "dropped_web" and e.side == "lost" for e in diff.events)
    lofted = set(diff.shared)
    assert lofted == set(shared_loft_names([r.name for r in root.runs()], [r.name for r in tip.runs()]))
    lost_webs = {e.name for e in diff.events if e.kind == "dropped_web"}
    assert lost_webs
    assert lost_webs.isdisjoint(lofted)


def test_web_feet_on_stacked_loop():
    cut = _geo().section_at(4.0)
    g = cut.shell_graph()
    hosted = g.with_web_feet_on_loop()
    loop_xy = np.vstack([r.points[:, :2] for r in g.loop])
    for web in hosted.webs:
        for foot in (web.points[0, :2], web.points[-1, :2]):
            assert float(np.min(np.linalg.norm(loop_xy - foot, axis=1))) == pytest.approx(0.0)


def test_shell_strips_from_graph_no_extra_cap_surface():
    from implicit_geometry_kernel import shell_strips_from_graph

    cut = _geo().section_at(4.0)
    g = cut.shell_graph()
    strips = shell_strips_from_graph(cut)
    assert len(strips) == len(g.runs())
    assert len(strips) < len(cut.walls)
    cap = next(s for s in strips if s.kind == "cap")
    run = next(r for r in g.loop if r.name == cap.label)
    assert cap.thickness_m == run.thickness
    assert cap.requires_eccentricity
    assert cap.eccentricity_m == run.eccentricity_m
    assert cap.midline_twisted_pitch_axis_blade_frame.shape[1] == 2
    assert any(w.kind == "cap" for w in cut.walls)


def test_rib_face_at_member_z_only():
    rib = TransverseMember("rib_mid", 5.0, 0.012)
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                        shear_web_thickness=0.002, web_positions=(0.18, 0.48),
                        transverse=(rib,), pivot_fraction=0.3)
    assert geo.forced_stations() == (5.0,)
    at_rib = geo.section_at(5.0)
    assert at_rib.transverse == (rib,)
    faces = at_rib.shell_graph().faces
    assert faces
    assert all(f.surface == "rib" and f.kind == "rib_mid" for f in faces)
    assert all(np.allclose(f.points[0], f.points[-1]) for f in faces)
    assert all(f.points.shape[1] == 3 and np.allclose(f.points[:, 2], 5.0) for f in faces)
    oml = geo.oml_sdf()
    for face in faces:
        phi = oml.evaluate(face.points).phi
        assert np.all(phi < 1e-6)
    other = geo.section_at(2.0).shell_graph()
    assert other.faces == ()
    assert all(r.surface != "rib" for r in other.runs())
    assert all(r.surface != "rib" for r in at_rib.shell_graph().runs())
    diff = topology_diff(at_rib.shell_graph(), other)
    assert any(e.kind == "lost_rib" and e.side == "lost" and e.name.startswith("rib_mid")
               for e in diff.events)
    gained = topology_diff(other, at_rib.shell_graph())
    assert any(e.kind == "gained_rib" and e.side == "gained" for e in gained.events)
    assert all(not n.startswith("rib_mid") for n in diff.shared)


def test_partial_span_stringer_branch_and_topology():
    stringer = Stringer("str_le", position=0.12, thickness=0.008, depth=0.04,
                        extent=(2.0, 6.0), side="pressure")
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                        shear_web_thickness=0.002, web_positions=(0.18, 0.48),
                        stringers=(stringer,), pivot_fraction=0.3)
    assert geo.forced_stations() == (2.0, 6.0)
    inside = geo.section_at(4.0).shell_graph()
    assert len(inside.branches) == 1
    run = inside.branches[0]
    assert run.surface == "stringer" and run.name == "str_le"
    assert run.points.shape == (2, 3)
    assert float(np.linalg.norm(run.points[1, :2] - run.points[0, :2])) == pytest.approx(0.04, rel=0.2)
    assert geo.section_at(0.5).shell_graph().branches == ()
    assert geo.section_at(8.0).shell_graph().branches == ()
    assert all(r.surface != "stringer" for r in inside.runs())
    diff = topology_diff(inside, geo.section_at(8.0).shell_graph())
    assert any(e.kind == "lost_stringer" and e.name == "str_le" for e in diff.events)
    assert "str_le" not in diff.shared


@pytest.mark.parametrize("key", SYSTEM_TYPE_LAYOUT_KEYS)
def test_wall_mid_matches_graph_run_interiors(key):
    """Exactness: ``wall_mid_sdf`` is a view of ``shell_graph`` (all 40 layouts)."""
    stations = [
        Station(0.0, 1.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(3.0, 0.8, Airfoil.naca4("2418"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations,
        key,
        skin_thickness=0.004,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005, pivot_fraction=0.3)
    try:
        cut = geo.section_at(1.5)
    except InfeasibleSectionGeometry:
        pytest.skip("layout infeasible at these gauges")
    mid = geo.wall_mid_sdf()
    worst = 0.0
    for run in cut.shell_graph().loop:
        if len(run.points) < 2:
            continue
        p = run.points[len(run.points) // 2].reshape(1, 3)
        worst = max(worst, abs(float(mid.evaluate(p).phi[0])))
    assert worst < 5e-4, f"{key} worst |phi| at loop-run interiors: {worst}"
