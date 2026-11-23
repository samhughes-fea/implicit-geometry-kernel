"""Sub-partitions step properties, never geometry (decision 0016)."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.mesh import shell_mesh
from implicit_geometry_kernel.mesh.shell.fields import assign_sub_partitions, element_tz
from implicit_geometry_kernel.structure import GridThickness, SubPartition, TableLaw, check_sub_partitions


def _blade():
    foil = Airfoil.naca4("0018")
    stations = [
        Station(0.0, 1.0, foil, twist=5.0, pitch_axis=0.3),
        Station(4.0, 0.7, foil, twist=0.0, pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-0R", 0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3,
    )


@pytest.fixture(scope="module")
def blade_and_mesh():
    geo = _blade()
    return geo, shell_mesh(geo, h=0.08)


def _skin_host(mesh):
    names = [p for p, s in zip(mesh.partition, mesh.surface) if s == "skin"]
    return max(set(names), key=names.count)


def test_sub_partition_validates_windows():
    with pytest.raises(ValueError):
        SubPartition("a", "skin", (0.6, 0.4), 0.01)
    with pytest.raises(ValueError):
        SubPartition("a", "skin", (0.1, 0.4), 0.01, z_range=(2.0, 1.0))
    with pytest.raises(ValueError):
        SubPartition("a", "skin", (0.1, 0.4), -0.01)
    s = SubPartition("a", "skin", (0.1, 0.4), 0.01, z_range=(1.0, 2.0))
    assert s.contains([0.1, 0.4, 0.5], 1.5).tolist() == [True, True, False]
    assert not bool(s.contains(0.2, 2.5))


def test_thickness_forms():
    z = [0.0, 4.0]
    assert SubPartition("a", "h", (0, 1), 0.01).thickness_at(0.3, 2.0) == pytest.approx(0.01)
    assert SubPartition("a", "h", (0, 1), TableLaw(z, [0.01, 0.03])).thickness_at(0.3, 2.0) == pytest.approx(0.02)
    grid = GridThickness([0.0, 1.0], z, [[0.01, 0.01], [0.03, 0.03]])
    assert SubPartition("a", "h", (0, 1), grid).thickness_at(0.5, 2.0) == pytest.approx(0.02)


def test_overlap_and_duplicate_checks():
    a = SubPartition("a", "h", (0.0, 0.5), 0.01)
    b = SubPartition("b", "h", (0.5, 1.0), 0.01)
    assert check_sub_partitions([a, b]) == []
    c = SubPartition("c", "h", (0.4, 0.6), 0.01)
    assert any("overlap" in m for m in check_sub_partitions([a, c]))
    d = SubPartition("d", "h", (0.4, 0.6), 0.01, z_range=(5.0, 6.0))
    e = SubPartition("e", "h", (0.4, 0.6), 0.01, z_range=(1.0, 2.0))
    assert check_sub_partitions([d, e]) == []
    assert any("duplicate" in m for m in check_sub_partitions([a, a]))


def test_assign_steps_properties_not_nodes(blade_and_mesh):
    geo, mesh = blade_and_mesh
    host = _skin_host(mesh)
    nodes_before = np.array(mesh.nodes, float, copy=True)
    band = SubPartition("flange", host, (0.2, 0.5), 0.02, z_range=(1.0, 3.0),
                        eccentricity=(0.005, 0.001), kind="bond_flange")
    out = assign_sub_partitions(mesh, geo, [band])
    np.testing.assert_array_equal(np.asarray(mesh.nodes, float), nodes_before)
    t, z = element_tz(mesh, geo)
    on = np.array([p == host for p in mesh.partition])
    inside = on & (t >= 0.2) & (t <= 0.5) & (z >= 1.0) & (z <= 3.0)
    assert inside.any()
    assert (out.band[inside] == "flange").all()
    assert (out.band[~inside] == "").all()
    assert out.thickness[inside] == pytest.approx(0.02)
    np.testing.assert_allclose(out.thickness[~inside], np.asarray(mesh.thickness)[~inside])
    np.testing.assert_allclose(out.eccentricity[inside], [[0.005, 0.001]] * int(inside.sum()))
    np.testing.assert_allclose(out.eccentricity[~inside, 0], np.asarray(mesh.eccentricity)[~inside])
    assert out.straddles.any()
    assert not out.straddles[~on].any()


def test_assign_rejects_overlap(blade_and_mesh):
    geo, mesh = blade_and_mesh
    host = _skin_host(mesh)
    subs = [SubPartition("a", host, (0.0, 0.5), 0.01), SubPartition("b", host, (0.4, 0.6), 0.01)]
    with pytest.raises(ValueError, match="overlap"):
        assign_sub_partitions(mesh, geo, subs)


def test_assign_with_no_bands_is_identity(blade_and_mesh):
    geo, mesh = blade_and_mesh
    out = assign_sub_partitions(mesh, geo, [])
    np.testing.assert_allclose(out.thickness, mesh.thickness)
    assert (out.band == "").all()
    assert not out.straddles.any()


def test_caps_are_whole_run_sub_partitions(blade_and_mesh):
    geo, mesh = blade_and_mesh
    graph = geo.section_at(2.0).shell_graph()
    caps = graph.cap_sub_partitions()
    cap_runs = [r for r in graph.loop if r.kind == "cap"]
    assert len(caps) == len(cap_runs) > 0
    for sub, run in zip(caps, cap_runs):
        assert sub.host == run.name and sub.kind == "spar_cap"
        assert sub.t_range == (0.0, 1.0)
        assert sub.thickness_at(0.5, 2.0) == pytest.approx(run.thickness)
        assert sub.eccentricity == (pytest.approx(run.eccentricity_m), 0.0)
    out = assign_sub_partitions(mesh, geo, caps)
    on_cap = np.array([p in {c.host for c in caps} for p in mesh.partition])
    np.testing.assert_allclose(out.thickness, mesh.thickness)
    np.testing.assert_allclose(out.eccentricity[:, 0], mesh.eccentricity)
    assert (out.band[on_cap] != "").all()
    assert not out.straddles.any()
