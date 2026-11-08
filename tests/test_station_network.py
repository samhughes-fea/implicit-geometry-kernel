"""Station network and wall-meta field from section_at."""

import numpy as np

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, station_network, wall_meta_from_cut
from implicit_geometry_kernel.sdf import FeatureTable, GeometrySolver


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                         shear_web_thickness=0.002, web_positions=(0.18, 0.48), pivot_fraction=0.3)


def test_station_network_ids_stable_across_span():
    geo = _geo()
    table = FeatureTable()
    n0 = station_network(geo.section_at(2.0), table)
    n1 = station_network(geo.section_at(8.0), table)
    assert n0.loop_features == n1.loop_features
    assert {w.kind for w in n0.walls} <= {"skin", "web", "cap"}
    assert {"skin", "web", "cap"} <= {w.kind for w in n0.walls}
    assert [table.by_name(w.label).id for w in n0.walls] == [w.feature_id for w in n0.walls]
    assert len(table) == len(n0.walls)


def test_wall_meta_tags_a_web_midpoint():
    geo = _geo()
    cut = geo.section_at(4.0)
    field, table = wall_meta_from_cut(cut)
    web = next(w for w in cut.walls if w.kind == "web")
    mid = web.xy[len(web.xy) // 2]
    sample = field.evaluate_plane(mid.reshape(1, 2), cut.z)
    assert sample.phi[0] < 0.0
    assert sample.meta[0] == table.by_name(web.label).id


def test_oml_sdf_project_stays_in_plane():
    geo = _geo()
    field = geo.oml_sdf(meta=0)
    ring = geo.oml_polyline(4.0)
    outside = ring[0] + np.array([0.05, 0.0, 0.0])
    out = GeometrySolver(tolerance=1e-7, max_step=0.05).project(field, outside.reshape(1, 3))
    assert out.converged[0]
    np.testing.assert_allclose(out.points[0, 2], ring[0, 2], atol=1e-12)
    assert abs(field.evaluate(out.points).phi[0]) < 1e-6
