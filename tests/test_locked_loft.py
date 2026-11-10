"""Zone-locked loft: signature bands, shared partition vertices, write_mid_surfaces."""

from pathlib import Path

import numpy as np

from implicit_geometry_kernel import (
    Airfoil,
    BladeGeometry,
    LockedLoftPlan,
    Station,
    resample_locked,
    shell_loft_zones,
)
from implicit_geometry_kernel.blade.loft.locked_loft import locked_mid_surface_patches
from implicit_geometry_kernel.blade.shell_graph import partition_edge_curve


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                         shear_web_thickness=0.002, web_positions=(0.18, 0.48), pivot_fraction=0.3)


def test_shell_loft_zones_break_on_dropped_web():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, web_offsets=(0.4, 1.2), pivot_fraction=0.3)
    zs = np.linspace(0.0, 10.0, 21)
    graphs = [geo.section_at(float(z)).shell_graph() for z in zs]
    zones = shell_loft_zones(graphs)
    assert len(zones) >= 2
    sigs = [tuple(r.name for r in graphs[sl][0].webs) for sl in zones]
    assert len(set(sigs)) == len(sigs)


def test_shell_loft_zones_constant_signature_is_one_band():
    geo = _geo()
    graphs = [geo.section_at(z).shell_graph() for z in (2.0, 4.0, 8.0)]
    zones = shell_loft_zones(graphs)
    assert zones == [slice(0, 3)]


def test_resample_locked_shares_partition_endpoints():
    g = _geo().section_at(4.0).shell_graph()
    plan = LockedLoftPlan.from_graph(g, n_chord=12)
    samples = resample_locked(g, plan)
    assert plan.names
    edges = g.partition_edges()
    assert edges
    for edge in edges:
        pa, pb = samples[edge.name_a], samples[edge.name_b]
        da = min(float(np.linalg.norm(pa[0] - edge.point)),
                 float(np.linalg.norm(pa[-1] - edge.point)))
        db = min(float(np.linalg.norm(pb[0] - edge.point)),
                 float(np.linalg.norm(pb[-1] - edge.point)))
        assert da == 0.0 and db == 0.0
        ends_a = (pa[0], pa[-1])
        ends_b = (pb[0], pb[-1])
        assert any(np.allclose(ea, eb) for ea in ends_a for eb in ends_b)


def test_resample_locked_web_feet_on_loop_samples():
    g = _geo().section_at(4.0).shell_graph()
    plan = LockedLoftPlan.from_graph(g, n_chord=12)
    samples = resample_locked(g, plan)
    loop = np.vstack([samples[n] for n in plan.loop_names])
    for name in plan.web_names:
        for i in (0, -1):
            d = float(np.min(np.linalg.norm(loop - samples[name][i], axis=1)))
            assert d == 0.0


def test_moving_cap_windows_joint_matches_partition_edge_curve():
    def windows(z):
        return ((0.20 + 0.01 * float(z), 0.50 + 0.01 * float(z)),)

    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
        shear_web_thickness=0.002, web_positions=(0.18, 0.48),
        cap_windows=windows, pivot_fraction=0.3)
    zs = [2.0, 4.0, 8.0]
    graphs = [geo.section_at(z).shell_graph() for z in zs]
    assert len(shell_loft_zones(graphs)) == 1
    plan = LockedLoftPlan.from_graph(graphs[0], n_chord=16)
    cap = next(e for e in graphs[0].partition_edges()
               if "spar_cap" in e.name_a or "spar_cap" in e.name_b)
    curve = partition_edge_curve(geo, cap.name_a, cap.name_b, zs)
    assert len(curve) == 3
    for g, xyz in zip(graphs, curve):
        samples = resample_locked(g, plan)
        pa, pb = samples[cap.name_a], samples[cap.name_b]
        joint = pa[_closest_end(pa, xyz)]
        np.testing.assert_allclose(joint, xyz)
        np.testing.assert_allclose(pb[_closest_end(pb, xyz)], xyz)


def _closest_end(pts, point):
    d0 = float(np.linalg.norm(pts[0] - point))
    d1 = float(np.linalg.norm(pts[-1] - point))
    return 0 if d0 <= d1 else -1


def test_write_mid_surfaces_obj_shares_seam_vertices(tmp_path: Path):
    geo = _geo()
    zs = np.linspace(geo.stations[0].z, geo.stations[-1].z, 5)
    graphs = [geo.shell_graph(float(z)) for z in zs]
    patches = {name: rings for name, rings, _surface in locked_mid_surface_patches(graphs, n_chord=12)}
    g = graphs[0]
    edge = next(e for e in g.partition_edges() if "spar_cap" in e.name_a or "spar_cap" in e.name_b)
    a, b = patches[edge.name_a], patches[edge.name_b]
    assert a.shape[0] == b.shape[0]
    for ring_a, ring_b in zip(a, b):
        ends_a = np.vstack([ring_a[0], ring_a[-1]])
        ends_b = np.vstack([ring_b[0], ring_b[-1]])
        dist = float(np.linalg.norm(ends_a[:, None, :] - ends_b[None, :, :], axis=2).min())
        assert dist == 0.0, dist
    obj = Path(geo.write_mid_surfaces(tmp_path / "mids.obj", n_span=5, n_chord=12)).read_text(encoding="utf-8")
    groups = [line[2:] for line in obj.splitlines() if line.startswith("g ")]
    assert groups
    assert any(n.startswith("cap_") or n.startswith("bay_") or n.startswith("te") for n in groups)


def test_shared_names_loft_through_web_drop():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, web_offsets=(0.4, 1.2), pivot_fraction=0.3)
    zs = np.linspace(0.0, 10.0, 11)
    graphs = [geo.section_at(float(z)).shell_graph() for z in zs]
    assert len(shell_loft_zones(graphs)) >= 2
    patches = locked_mid_surface_patches(graphs, n_chord=8)
    first = {r.name for r in graphs[0].loop}
    last = {r.name for r in graphs[-1].loop}
    shared = first & last
    assert shared
    by_name: dict[str, list] = {}
    for name, rings, _surface in patches:
        by_name.setdefault(name, []).append(rings)
    for name in shared:
        present = sum(1 for g in graphs if any(r.name == name for r in g.runs()))
        covered = sum(r.shape[0] for r in by_name[name])
        assert covered == present, (name, covered, present)


def test_locked_patches_split_when_web_drops():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, web_offsets=(0.4, 1.2), pivot_fraction=0.3)
    zs = np.linspace(0.0, 10.0, 11)
    graphs = [geo.section_at(float(z)).shell_graph() for z in zs]
    patches = locked_mid_surface_patches(graphs, n_chord=8)
    web_names = {p[0] for p in patches if "web" in p[0]}
    assert web_names
    spans = {p[0]: p[1].shape[0] for p in patches}
    assert any(n < len(graphs) for n in spans.values())
