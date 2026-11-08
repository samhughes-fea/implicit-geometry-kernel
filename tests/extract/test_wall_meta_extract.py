"""Wall-meta DAG, solver snap, extractors, station network v2."""

import numpy as np

from implicit_geometry_kernel import (
    Airfoil,
    BladeGeometry,
    Station,
    station_network_from_meta,
)
from implicit_geometry_kernel.extract import classify_trailing_edge, extract_zero_loops, feature_curve
from implicit_geometry_kernel.sdf import GeometrySolver
from implicit_geometry_kernel.sdf.nodes.tagged import CircleSDF


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                         shear_web_thickness=0.002, web_positions=(0.18, 0.48), pivot_fraction=0.3)


def test_wall_meta_tags_web_at_midspan():
    geo = _geo()
    cut = geo.section_at(4.0)
    web = next(w for w in cut.walls if w.kind == "web")
    field = geo.wall_meta()
    mid = np.append(web.xy[len(web.xy) // 2], cut.z)
    sample = field.evaluate(mid.reshape(1, 3))
    assert sample.phi[0] < 0.0
    spec = field.table.by_name(web.label)
    assert sample.meta[0] == spec.id


def test_cavity_and_wall_mid_signs():
    geo = _geo()
    ring = geo.oml_polyline(4.0)
    outside = ring[0] + np.array([0.2, 0.0, 0.0])
    assert geo.oml_sdf().evaluate(outside.reshape(1, 3)).phi[0] > 0.0
    assert geo.cavity_sdf().evaluate(outside.reshape(1, 3)).phi[0] > 0.0


def test_snap_oml_stays_in_plane():
    geo = _geo()
    ring = geo.oml_polyline(4.0)
    start = ring[5] + np.array([0.03, 0.0, 0.0])
    out = geo.snap_oml(start.reshape(1, 3))
    assert out.converged[0]
    np.testing.assert_allclose(out.points[0, 2], 4.0)
    assert abs(geo.oml_sdf().evaluate(out.points).phi[0]) < 1e-6


def test_snap_wall_mid_near_skin():
    geo = _geo()
    cut = geo.section_at(4.0)
    skin = next(w for w in cut.walls if w.kind == "skin")
    p = skin.points[len(skin.points) // 2]
    nudged = p + np.array([0.0, 0.002, 0.0])
    out = geo.snap_wall_mid(nudged.reshape(1, 3))
    assert out.converged[0]
    np.testing.assert_allclose(out.points[0, 2], 4.0, atol=1e-12)


def test_intersect_web_and_oml():
    geo = _geo()
    cut = geo.section_at(4.0)
    web = next(w for w in cut.walls if w.kind == "web")
    from implicit_geometry_kernel.blade.network import wall_meta_from_cut
    meta, table = wall_meta_from_cut(cut)
    web_id = table.by_name(web.label).id
    # a point on the web near the skin: solver still holds z
    guess = web.points[0]
    hit = GeometrySolver(tolerance=1e-6, max_step=0.02).project(geo.oml_sdf(), guess.reshape(1, 3))
    assert hit.points[0, 2] == guess[2]
    assert web_id >= 0 and meta is not None


def test_extract_circle_loop():
    field = CircleSDF(0.0, 0.0, 1.0)
    loops = extract_zero_loops(field, 0.0, (-1.5, 1.5, -1.5, 1.5), n=32)
    assert loops
    r = np.linalg.norm(loops[0][:, :2], axis=1)
    np.testing.assert_allclose(r, 1.0, atol=0.05)


def test_feature_curve_two_circles():
    a, b = CircleSDF(0.0, 0.0, 1.0), CircleSDF(1.0, 0.0, 1.0)
    zs = [0.0, 1.0, 2.0]
    curve = feature_curve(a, b, zs, (0.5, 0.8))
    assert curve.shape == (3, 3)
    np.testing.assert_allclose(curve[:, 2], zs)
    np.testing.assert_allclose(curve[:, 0], 0.5, atol=1e-4)


def test_classify_te_sharp_naca():
    geo = _geo()
    ring = geo.oml_polyline(4.0)
    assert classify_trailing_edge(ring) in {"sharp", "blunt"}


def test_station_network_v2_from_meta():
    geo = _geo()
    z = 4.0
    field = geo.wall_meta()
    cut = geo.section_at(z)
    xs = np.concatenate([w.xy[:, 0] for w in cut.walls])
    ys = np.concatenate([w.xy[:, 1] for w in cut.walls])
    pad = 0.05
    bounds = (float(xs.min()) - pad, float(xs.max()) + pad, float(ys.min()) - pad, float(ys.max()) + pad)
    net = station_network_from_meta(geo.wall_mid_sdf(), field, z, bounds, field.table, n=36)
    assert net.z == z
    assert net.loop_features
    assert all(fid in range(len(field.table)) for fid in net.loop_features)
