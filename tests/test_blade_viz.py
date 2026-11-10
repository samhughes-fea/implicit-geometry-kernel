"""Example viz helpers: suction-side names and twist-aware web fractions."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples"))

from _shared.blade_viz import (  # noqa: E402
    _CSG_FIG_H_MIN,
    _CSG_FIG_W,
    _CSG_QT_DEPTH,
    _attach_union_phi,
    _blade_xy_to_host_naca,
    _conformal_samples,
    _csg_mesh_figsize,
    _csg_quadtree_spec,
    _leaf_rects,
    _map_naca_chord_xy,
    _property_quadtree_spec,
    _sample_property_qt,
    _sample_wall_meta,
    _station_limits,
    _te_mesh_limits,
    constructive_polygons_blade_xy,
    is_suction_role,
    live_web_fracs,
    sample_planform,
)
from implicit_geometry_kernel.integration.quadtree import leaf_bounds  # noqa: E402
from implicit_geometry_kernel.section.boundary.labels import (  # noqa: E402
    is_core_label,
    is_wall_solid_label,
)
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.blade.topology import _chord_fracs  # noqa: E402


def test_is_suction_role_matches_authored_loop_names():
    assert is_suction_role("bay_su_0")
    assert is_suction_role("bay_su_1")
    assert is_suction_role("te_su")
    assert is_suction_role("le_su")
    assert is_suction_role("cap_su")
    assert is_suction_role("spar_cap_suction_surface")
    assert not is_suction_role("bay_pr_0")
    assert not is_suction_role("te_pr")
    assert not is_suction_role("fore_web")
    assert not is_suction_role("web_0")
    assert not is_suction_role("tip_closure:cell_0")
    assert not is_suction_role("rib_0")


def test_live_web_fracs_match_templates_along():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), twist=12.0, pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), twist=8.0, pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    cut = geo.section_at(5.0)
    live = live_web_fracs(cut)
    expected = []
    for run in cut.shell_graph().webs:
        pts = np.asarray(run.points, float)[:, :2]
        expected.append(float(np.clip(_chord_fracs(cut, pts[:1])[0], 0.0, 1.0)))
    assert live == expected
    assert live
    naive = []
    for run in cut.shell_graph().webs:
        xy = np.asarray(run.points, float)[:, :2]
        naive.append(float(np.mean(xy[:, 0] / max(cut.chord, 1e-9))))
    assert naive != live


def test_sample_planform_includes_rib_lines():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-4R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    data = sample_planform(geo, n=11)
    assert len(data["ribs"]) == 4
    zs = [row[0] for row in data["ribs"]]
    assert min(zs) > 0.0
    assert max(zs) == pytest.approx(10.0)
    for z, x_le, x_te, name in data["ribs"]:
        assert x_te > x_le
        assert name


def test_csg_wall_meta_quadtree_caps_leaf_depth_at_ten():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), twist=12.0, pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), twist=8.0, pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3,
    )
    cut = geo.section_at(5.0)
    xlim, ylim = _station_limits(cut)
    bounds, mesh = _sample_wall_meta(geo, cut, cut.z, xlim, ylim)
    levels = np.asarray(mesh._leaf_levels, int)
    chord = float(cut.chord)
    t = min(float(w.thickness) for w in cut.walls if float(w.thickness) > 0.0)
    radius = 0.5 * t
    assert bounds.shape[0] >= 1
    assert int(levels.max()) <= _CSG_QT_DEPTH
    assert int(mesh.spec.max_depth) == 10
    assert _CSG_QT_DEPTH == 10
    assert mesh.spec.thin_band is True
    assert mesh.spec.min_cell_size == pytest.approx(min(0.0015 * chord, 0.5 * radius))
    assert mesh.spec.refine_band == pytest.approx(radius)
    phi = np.asarray(mesh.phi_by_label["wall"], float)
    assert np.any(phi <= 0.0)
    width = max(float(xlim[1] - xlim[0]), 1e-9)
    height = max(float(ylim[1] - ylim[0]), 1e-9)
    dx = max(t / 4.0, 1e-6)
    n_uniform = (int(np.ceil(width / dx)) + 1) * (int(np.ceil(height / dx)) + 1)
    assert bounds.shape[0] < n_uniform


def test_csg_mesh_figsize_follows_section_aspect():
    wide, tall = _csg_mesh_figsize((-1.0, 1.0), (-1.0, 1.0))
    assert wide == pytest.approx(_CSG_FIG_W)
    assert tall == pytest.approx(2.0 * (_CSG_FIG_W / 3.0) + 2.2)
    thin_w, thin_h = _csg_mesh_figsize((0.0, 1.0), (-0.05, 0.05))
    assert thin_w == pytest.approx(_CSG_FIG_W)
    assert thin_h < 10.0
    assert thin_h >= _CSG_FIG_H_MIN


def _foil_cut():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), twist=12.0, pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), twist=8.0, pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3,
    )
    return geo, geo.section_at(5.0)


def _shoelace(xy: np.ndarray) -> float:
    p = np.asarray(xy, float)[:, :2]
    if len(p) < 3:
        return 0.0
    if not np.allclose(p[0], p[-1]):
        p = np.vstack([p, p[0]])
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))))


def test_csg_constructive_polygons_are_wall_solids_not_oml_disk():
    geo, cut = _foil_cut()
    polys = constructive_polygons_blade_xy(cut)
    keys = set(polys)
    assert any(is_core_label(k) for k in keys)
    assert any("skin" in k for k in keys)
    assert any("shear_web" in k or k.startswith("web") for k in keys)
    assert any("spar_cap" in k for k in keys)
    wall_area = 0.0
    for key, rings in polys.items():
        if not is_wall_solid_label(key):
            continue
        wall_area += _shoelace(rings[0])
        for hole in rings[1:]:
            wall_area -= _shoelace(hole)
    oml = np.asarray(geo.oml_polyline(float(cut.z)), float)
    oml_area = _shoelace(oml)
    assert wall_area > 0.0
    assert wall_area < 0.45 * oml_area
    skin = next(v for k, v in polys.items() if "skin_suction" in k)
    sx = np.vstack(skin)
    assert sx[:, 0].min() == pytest.approx(float(oml[:, 0].min()), abs=0.02 * cut.chord)
    assert sx[:, 0].max() == pytest.approx(float(oml[:, 0].max()), abs=0.02 * cut.chord)


def test_csg_property_quadtree_uses_thin_band_wall_tubes():
    geo, cut = _foil_cut()
    assert abs(float(cut.twist)) > 1.0
    spec = _property_quadtree_spec(geo, cut, cut.z)
    occ = _csg_quadtree_spec(geo, cut, cut.z)
    assert spec.thin_band is True
    assert occ.thin_band is True
    assert spec.max_depth == occ.max_depth == _CSG_QT_DEPTH
    assert spec.min_cell_size == pytest.approx(occ.min_cell_size)
    assert spec.refine_band == pytest.approx(occ.refine_band)
    k = float(getattr(cut, "scale", 1.0) or 1.0)
    naca = np.array([[0.5 * float(cut.chord) * k, 0.0]], float)
    blade = _map_naca_chord_xy(cut, naca)
    back = _blade_xy_to_host_naca(cut, blade)
    assert back == pytest.approx(naca, abs=1e-9)
    mesh = _sample_property_qt(geo, cut, cut.z)
    assert mesh.spec.thin_band is True
    assert int(mesh.spec.max_depth) == 10
    rects = _leaf_rects(leaf_bounds(mesh))
    assert len(rects) >= 1
    for quad in rects[:32]:
        xs, ys = quad[:, 0], quad[:, 1]
        assert len(np.unique(np.round(xs, 12))) == 2
        assert len(np.unique(np.round(ys, 12))) == 2
    assert _attach_union_phi(mesh)
    phi = np.asarray(mesh.phi_by_label["_union_walls"], float)
    assert np.any(phi <= 0.0)
    inside = float((phi[:, :4] <= 0.0).mean())
    assert 0.0 < inside < 0.5
    labels = {str(k) for k in mesh.phi_by_label if k != "_union_walls"}
    assert any(is_core_label(k) for k in labels)
    assert any(is_wall_solid_label(k) for k in labels)


def test_te_mesh_limits_contain_te_and_are_smaller_than_station():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.004,
        pivot_fraction=0.3, spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )
    cut = geo.section_at(5.0)
    samples = _conformal_samples(geo, cut, 5.0)
    te_x = []
    for name in ("te_su", "te_pr"):
        if name in samples and len(samples[name]):
            te_x.append(float(np.asarray(samples[name], float)[-1, 0]))
    assert te_x
    xlim, ylim = _te_mesh_limits(cut, samples)
    full_x, full_y = _station_limits(cut)
    for x in te_x:
        assert xlim[0] < x < xlim[1]
    assert (xlim[1] - xlim[0]) < (full_x[1] - full_x[0])
    assert (ylim[1] - ylim[0]) == pytest.approx(xlim[1] - xlim[0])
