"""Thickness as a field t = f(s, z) on each partition (decision 0014), properties only (0016)."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.mesh import shell_mesh
from implicit_geometry_kernel.mesh.shell.fields import element_tz, node_thickness, sample_thickness
from implicit_geometry_kernel.structure import (
    GridThickness,
    TableLaw,
    ThicknessField,
    UnderSampledSpanLaw,
)


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
    # Ruled: its leading-edge counts follow curvature, so h=0.08 resolves the nose.
    # The partition mesher spaces t evenly and needs about h=0.04 for the same check.
    return geo, shell_mesh(geo, h=0.08, method="ruled")


def test_grid_thickness_is_linear_across_and_along():
    field = GridThickness([0.0, 1.0], [0.0, 4.0], [[0.01, 0.02], [0.03, 0.04]])
    assert field(0.5, 2.0) == pytest.approx(0.025)
    assert field([0.0, 1.0], 0.0).tolist() == pytest.approx([0.01, 0.03])
    assert field(2.0, 9.0) == pytest.approx(0.04)


def test_grid_thickness_pchip_rows_are_sampling_checked():
    z = np.linspace(0.0, 4.0, 4)
    with pytest.raises(UnderSampledSpanLaw):
        GridThickness([0.0, 1.0], z, np.vstack([0.01 + 0.005 * np.sin(z), 0.01 + 0.0 * z]), z_interp="pchip")
    zz = np.linspace(0.0, 4.0, 41)
    field = GridThickness([0.0, 1.0], zz, np.vstack([0.01 + 0.005 * np.sin(zz), 0.01 + 0.0 * zz]),
                          z_interp="pchip")
    assert field(0.0, 1.3) == pytest.approx(0.01 + 0.005 * np.sin(1.3), abs=1e-6)


def test_grid_thickness_rejects_bad_input():
    with pytest.raises(ValueError):
        GridThickness([0.0, 1.0], [0.0, 1.0], [[0.01, 0.02]])
    with pytest.raises(ValueError):
        GridThickness([0.0, 1.0], [0.0, 1.0], [[0.01, -0.02], [0.0, 0.0]])


def test_span_field_ignores_t():
    field = ThicknessField.from_span(TableLaw([0.0, 4.0], [0.01, 0.02]))
    assert field([0.0, 0.7], 2.0).tolist() == pytest.approx([0.015, 0.015])


def test_element_t_lies_in_the_unit_interval_on_skin_and_webs(blade_and_mesh):
    geo, mesh = blade_and_mesh
    t, z = element_tz(mesh, geo)
    on_patch = np.array([s in ("skin", "web") for s in mesh.surface])
    assert np.all(np.isfinite(t[on_patch]))
    assert np.all((t[on_patch] >= 0.0) & (t[on_patch] <= 1.0))
    assert np.all(np.isnan(t[~on_patch]))
    cap = [i for i, p in enumerate(mesh.partition) if p == "spar_cap_suction_surface"]
    assert t[cap].max() - t[cap].min() > 0.8


def test_fields_set_thickness_on_their_partition_only(blade_and_mesh):
    geo, mesh = blade_and_mesh
    field = GridThickness([0.0, 1.0], [0.0, 4.0], [[0.008, 0.004], [0.012, 0.006]])
    th = sample_thickness(mesh, geo, {"spar_cap_suction_surface": field})
    cap = np.array([p == "spar_cap_suction_surface" for p in mesh.partition])
    assert np.allclose(th[~cap], np.asarray(mesh.thickness)[~cap])
    t, z = element_tz(mesh, geo)
    assert np.allclose(th[cap], field(t[cap], z[cap]))
    assert np.array_equal(sample_thickness(mesh, geo, {}), np.asarray(mesh.thickness))


def test_node_thickness_does_not_average_across_partitions(blade_and_mesh):
    _, mesh = blade_and_mesh
    everywhere = node_thickness(mesh)
    assert np.all(np.isfinite(everywhere))
    cap_side = node_thickness(mesh, partition="spar_cap_suction_surface")
    cap_value = float(np.asarray(mesh.thickness)[list(mesh.partition).index("spar_cap_suction_surface")])
    on_cap = np.isfinite(cap_side)
    assert on_cap.sum() > 0
    assert np.allclose(cap_side[on_cap], cap_value)
    assert np.all(np.isnan(cap_side[~on_cap]))


def test_element_curvature_tags_every_element(blade_and_mesh):
    from implicit_geometry_kernel.mesh.shell.fields import element_curvature

    geo, mesh = blade_and_mesh
    curv = element_curvature(mesh, geo)
    assert "UNSET" not in set(curv["kind"])
    flat_planes = [i for i, s in enumerate(mesh.surface) if s in ("rib", "tip_closure")]
    assert all(curv["kind"][i] == "FLAT" and curv["k1"][i] == 0.0 for i in flat_planes)
    skin = [i for i, s in enumerate(mesh.surface) if s == "skin"]
    assert np.all(curv["k1"][skin] <= curv["k2"][skin])
    assert any(curv["kind"][i] != "FLAT" for i in skin)


def test_node_curvature_is_per_partition_and_matches_elements(blade_and_mesh):
    from implicit_geometry_kernel.mesh.shell.fields import element_curvature, node_curvature

    geo, mesh = blade_and_mesh
    part = "skin_suction_surface:leading_edge"
    nodes = node_curvature(mesh, geo, partition=part)
    elems = element_curvature(mesh, geo)
    members = [i for i, p in enumerate(mesh.partition) if p == part]
    on = sorted({k for i in members for k in mesh.elements[i]})
    off = sorted(set(range(len(mesh.nodes))) - set(on))
    assert np.all(np.isfinite(nodes["k2"][on]))
    assert np.all(np.isnan(nodes["k2"][off]))
    assert "UNSET" not in set(nodes["kind"][on])
    # The mean of an element's node curvatures is close to its centroid value.
    mean = np.array([np.mean(nodes["k2"][list(mesh.elements[i])]) for i in members])
    rel = np.abs(mean - elems["k2"][members]) / np.maximum(np.abs(elems["k2"][members]), 1e-9)
    assert float(np.median(rel)) < 0.1
    tip = [i for i, s in enumerate(mesh.surface) if s == "tip_closure"]
    flat = node_curvature(mesh, geo, partition=mesh.partition[tip[0]]) if tip else None
    if flat is not None and all(mesh.surface[i] == "tip_closure" for i, p in enumerate(mesh.partition) if p == mesh.partition[tip[0]]):
        assert set(flat["kind"][[k for k in mesh.elements[tip[0]]]]) == {"FLAT"}
    with pytest.raises(ValueError):
        node_curvature(mesh, geo, partition="nope")
