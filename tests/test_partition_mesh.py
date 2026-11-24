"""The (t, z) partition mesher beside the ruled-strip mesher (decision 0015)."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.mesh import shell_mesh
from implicit_geometry_kernel.mesh.shell.fields import assign_sub_partitions
from implicit_geometry_kernel.mesh.shell.partition_mesh import PartitionMeshError, partition_mesh
from implicit_geometry_kernel.structure import SubPartition

H = 0.08


def _blade(layout="2D-F-LC-0R", **kw):
    foil = Airfoil.naca4("0018")
    stations = [
        Station(0.0, 1.0, foil, twist=5.0, pitch_axis=0.3),
        Station(4.0, 0.7, foil, twist=0.0, pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, layout, 0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3, **kw,
    )


BAND_HOST = "skin_suction_surface:trailing_edge"


@pytest.fixture(scope="module")
def meshes():
    geo = _blade()
    band = SubPartition("flange", BAND_HOST, (0.13, 0.41), 0.006, kind="bond_flange")
    return geo, partition_mesh(geo, h=H, sub_partitions=[band]), shell_mesh(geo, h=H), band


def _area(mesh, surface):
    total = 0.0
    for e, s in zip(mesh.elements, mesh.surface):
        if s != surface:
            continue
        p = mesh.nodes[list(e)]
        if len(e) == 4:
            total += 0.5 * float(np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[1])))
        else:
            total += 0.5 * float(np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])))
    return total


def test_quads_on_skin_and_webs_triangles_on_the_tip(meshes):
    _geo, mesh, _ruled, _band = meshes
    assert set(mesh.surface) == {"skin", "web", "tip_closure"}
    for e, s, c in zip(mesh.elements, mesh.surface, mesh.constructor):
        if s in ("skin", "web"):
            assert len(e) == 4 and c == "partition"
        else:
            assert len(e) == 3 and c in ("inherited", "interior")
    tip = [i for i, s in enumerate(mesh.surface) if s == "tip_closure"]
    assert {mesh.partition[i] for i in tip} == {"skin"}
    assert all(mesh.z0[i] == max(mesh.stations) for i in tip)


def test_nodes_are_shared_by_identity_not_coordinate(meshes):
    _geo, mesh, _ruled, _band = meshes
    rounded = {tuple(np.round(p, 9)) for p in mesh.nodes}
    assert len(rounded) == len(mesh.nodes)
    used = {i for e in mesh.elements for i in e}
    assert used == set(range(len(mesh.nodes)))


def test_edges_are_conforming_with_junctions_at_web_feet(meshes):
    _geo, mesh, _ruled, _band = meshes
    use = mesh.edge_use()
    assert set(use.values()) <= {1, 2, 3}
    z0 = min(mesh.stations)
    for (a, b), n in use.items():
        za, zb = mesh.nodes[a, 2], mesh.nodes[b, 2]
        if n == 1:  # the tip is closed: open edges only at the root
            assert za == zb == z0
    spanwise_three = [(a, b) for (a, b), n in use.items() if n == 3 and mesh.nodes[a, 2] != mesh.nodes[b, 2]]
    assert len(spanwise_three) == 4 * (len(mesh.stations) - 1)  # one line per web foot, every bay


def test_skin_sheet_is_consistently_oriented(meshes):
    _geo, mesh, _ruled, _band = meshes
    directed: dict = {}
    for e, s in zip(mesh.elements, mesh.surface):
        if s != "skin":
            continue
        for a, b in zip(e, e[1:] + e[:1]):
            directed[(a, b)] = directed.get((a, b), 0) + 1
    assert max(directed.values()) == 1


def test_area_matches_the_ruled_strip_mesher(meshes):
    _geo, mesh, ruled, _band = meshes
    for surface in ("skin", "web", "tip_closure"):
        assert _area(mesh, surface) == pytest.approx(_area(ruled, surface), rel=1e-2)


def test_ribs_are_cells_on_the_station_nodes():
    geo = _blade(layout="2D-F-LC-2R")
    mesh = partition_mesh(geo, h=0.1)
    ruled = shell_mesh(geo, h=0.1)
    ribs = [i for i, s in enumerate(mesh.surface) if s == "rib"]
    assert {mesh.partition[i] for i in ribs} == {m.name for m in geo.transverse}
    for m in geo.transverse:
        assert any(abs(z - m.z) < 1e-9 for z in mesh.stations)
    # Coarse h: both meshers inscribe the rib outline in chords, the ruled one with more nodes.
    assert _area(mesh, "rib") == pytest.approx(_area(ruled, "rib"), rel=2e-2)
    # Each rib's boundary nodes are skin or web nodes of that station.
    for i in ribs:
        for k in mesh.elements[i]:
            if k < mesh.n_source:
                assert abs(mesh.nodes[k, 2] - mesh.z0[i]) < 1e-12
    assert len({tuple(np.round(p, 9)) for p in mesh.nodes}) == len(mesh.nodes)


def test_edges_across_partitions_are_at_most_h(meshes):
    _geo, mesh, _ruled, _band = meshes
    for e, s in zip(mesh.elements, mesh.surface):
        if s not in ("skin", "web"):
            continue
        p = mesh.nodes[list(e)]
        assert float(np.linalg.norm(p[1] - p[0])) <= H * 1.001 or float(np.linalg.norm(p[3] - p[0])) <= H * 1.001


def test_band_edges_are_mesh_lines(meshes):
    geo, mesh, _ruled, band = meshes
    out = assign_sub_partitions(mesh, geo, [band])
    assert (out.band == "flange").any()
    assert not out.straddles.any()
    assert out.thickness[out.band == "flange"] == pytest.approx(0.006)


def test_unknown_band_host_is_rejected():
    geo = _blade()
    with pytest.raises(ValueError, match="hosts"):
        partition_mesh(geo, h=0.2, zs=[0.0, 4.0], sub_partitions=[SubPartition("x", "nope", (0.1, 0.2), 0.01)])


def test_branches_are_refused_for_now():
    from implicit_geometry_kernel import Stringer

    geo = _blade(stringers=[Stringer("stiff", 0.6, 0.004, 0.02)])
    with pytest.raises(PartitionMeshError, match="branches"):
        partition_mesh(geo, h=0.2, zs=[0.0, 4.0])


def test_bad_size_is_rejected():
    with pytest.raises(ValueError):
        partition_mesh(_blade(), h=0.0)


def test_shell_mesh_method_switch():
    geo = _blade()
    mesh = shell_mesh(geo, h=0.2, zs=[0.0, 2.0, 4.0], method="partition")
    assert {c for c, s in zip(mesh.constructor, mesh.surface) if s in ("skin", "web")} == {"partition"}
    with pytest.raises(ValueError, match="needs h"):
        shell_mesh(geo, method="partition", zs=[0.0, 4.0])
    with pytest.raises(ValueError, match="method"):
        shell_mesh(geo, h=0.2, method="coons")
    auto = shell_mesh(geo, h=0.2, zs=[0.0, 2.0, 4.0], method="auto")
    assert {c for c, s in zip(auto.constructor, auto.surface) if s in ("skin", "web")} == {"partition"}
    with pytest.raises(ValueError, match="sub_partitions"):
        shell_mesh(geo, h=0.2, method="ruled", sub_partitions=[SubPartition("x", BAND_HOST, (0.1, 0.2), 0.01)])
    with pytest.raises(ValueError, match="sub_partitions"):
        shell_mesh(geo, h=0.2, h_over_c=0.1, sub_partitions=[SubPartition("x", BAND_HOST, (0.1, 0.2), 0.01)])


def test_auto_is_the_default_and_needs_sizing_by_h_alone():
    """Decision 0015 at 0.9.6: the partition mesher by default, the ruled mesher for its own controls."""
    geo = _blade()

    def sheet(mesh):
        return {c for c, s in zip(mesh.constructor, mesh.surface) if s in ("skin", "web")}

    default = shell_mesh(geo, h=0.2, zs=[0.0, 2.0, 4.0])
    assert sheet(default) == {"partition"}
    assert default.digest() == shell_mesh(geo, h=0.2, zs=[0.0, 2.0, 4.0], method="partition").digest()
    assert "partition" not in sheet(shell_mesh(geo, n_span=3))
    assert "partition" not in sheet(shell_mesh(geo, h=0.2, h_over_c=0.1))
    assert "partition" not in sheet(shell_mesh(geo, h=0.2, refine=[(1.0, 2.0, 0.1)]))


def test_auto_falls_back_to_the_ruled_mesh_where_the_partition_mesher_refuses():
    from implicit_geometry_kernel import Stringer

    geo = _blade(stringers=[Stringer("stiff", 0.6, 0.004, 0.02)])
    band = SubPartition("x", BAND_HOST, (0.1, 0.2), 0.01)
    mesh = shell_mesh(geo, h=0.2, zs=[0.0, 2.0, 4.0], method="auto", sub_partitions=[band])
    assert "partition" not in set(mesh.constructor)
    with pytest.raises(PartitionMeshError):
        shell_mesh(geo, h=0.2, zs=[0.0, 2.0, 4.0], method="partition")
