"""Two-view rule, loftable shell graph, merged wall_meta, stack-aware wall_mid."""

import numpy as np
import pytest

from implicit_geometry_kernel import (
    Airfoil,
    BladeGeometry,
    Station,
    TransverseMember,
    partition_edge_curve,
    shared_loft_names,
    station_network,
    topology_diff,
    wall_meta_from_cut,
)
from implicit_geometry_kernel.blade.cut import ShellGraph, ShellRun
from implicit_geometry_kernel.sdf import FeatureTable


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F", 0.003, spar_cap_thickness=0.008,
                         shear_web_thickness=0.002, web_positions=(0.18, 0.48))


def test_walls_keep_caps_mid_surface_does_not():
    cut = _geo().section_at(4.0)
    assert any(w.kind == "cap" for w in cut.walls)
    ms = cut.mid_surface()
    assert {s.kind for s in ms.segments} <= {"skin", "web"}
    assert "cap" not in {s.kind for s in ms.segments}


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
    t_cap = cut.wall(cap.name).thickness
    t_skin = cut.wall(host).thickness
    assert cap.thickness == t_cap + t_skin
    assert cap.layers == (cap.name, host)
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
    assert abs(cavity + oml) < 1e-3
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
    geo = BladeGeometry(stations, "2D-F", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
                        web_offsets=(0.4, 1.2))
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
    assert cap.midline_twisted_pitch_axis_blade_frame.shape[1] == 2
    assert any(w.kind == "cap" for w in cut.walls)


def test_rib_face_at_member_z_only():
    rib = TransverseMember("rib_mid", 5.0, 0.012)
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F", 0.003, spar_cap_thickness=0.008,
                        shear_web_thickness=0.002, web_positions=(0.18, 0.48),
                        transverse=(rib,))
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
