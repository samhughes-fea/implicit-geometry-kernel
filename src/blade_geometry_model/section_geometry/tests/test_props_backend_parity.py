"""Parity tests: grid vs quadtree section property backends."""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.api.factory import build_section_view, resolve_system_type
from blade_geometry_model.section_geometry.props_emit.area_mesh import build_shared_section_quadtree
from blade_geometry_model.section_geometry.props_emit.props_backend import quadtree_spec_for_airfoil
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.sdf.grid import SDFGrid
from blade_geometry_model.section_geometry.props_emit.export import SectionPropertiesReport
from blade_geometry_model.area_2d.quadtree import summarize_quadtree
from blade_geometry_model.section_geometry.sdf.cache import SectionEvalCache


@pytest.mark.parametrize("layout_key", ["2D-CN", "2D-F"])
def test_quadtree_backend_parity_with_grid(layout_key: str):
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type(layout_key), twist_angle_rad=0.0)
    grid = SDFGrid.from_airfoil(airfoil, nx=384, ny=180)
    cache = SectionEvalCache()
    grid_report = SectionPropertiesReport(
        section, grid, backend="grid", eval_cache=cache
    )
    qt_spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=8, max_cells=120000)
    qt_report = SectionPropertiesReport(
        section,
        grid,
        backend="quadtree",
        quadtree_spec=qt_spec,
        chord=1.0,
        twist_rad=0.0,
        airfoil_for_spec=airfoil,
    )
    for label in section.labels:
        g = grid_report[label]
        q = qt_report[label]
        if g["area"] < 1e-5:
            continue
        rtol = 0.12 if g["area"] < 0.002 else 0.08
        np.testing.assert_allclose(q["area"], g["area"], rtol=rtol, atol=5e-4)
        np.testing.assert_allclose(q["cx"], g["cx"], atol=0.05)
        np.testing.assert_allclose(q["cy"], g["cy"], atol=0.05)


def test_default_backend_is_hybrid():
    airfoil = AirfoilSDF.from_naca("0012", n_points=80, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    grid = SDFGrid.from_airfoil(airfoil, nx=200, ny=100)
    report = SectionPropertiesReport(section, grid)
    assert report.backend_meta.get("backend") == "hybrid"


def test_quadtree_root_cell_size_twist_invariant():
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    twist_rad = float(np.deg2rad(15.0))
    spec0 = quadtree_spec_for_airfoil(airfoil, chord=1.0)
    spec15 = quadtree_spec_for_airfoil(airfoil.rotate(twist_rad), chord=1.0)
    assert spec0.root_cell_size is not None
    assert spec15.root_cell_size is not None
    np.testing.assert_allclose(
        spec15.root_cell_size,
        spec0.root_cell_size,
        rtol=0.02,
        atol=1e-6,
    )


def test_quadtree_twist_resolution_matches():
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    layout = resolve_system_type("2D-CN")
    twist_rad = float(np.deg2rad(15.0))
    section0 = build_section_view(airfoil, layout, twist_angle_rad=0.0)
    section15 = build_section_view(airfoil, layout, twist_angle_rad=twist_rad)
    airfoil_b = airfoil.rotate(twist_rad)
    spec0 = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=7, max_cells=150000)
    spec15 = quadtree_spec_for_airfoil(airfoil_b, chord=1.0, max_depth=7, max_cells=150000)
    mesh0 = build_shared_section_quadtree(section0, spec0)
    mesh15 = build_shared_section_quadtree(section15, spec15, twist_rad=twist_rad)
    w0 = float(summarize_quadtree(mesh0)["min_cell_width"])
    w15 = float(summarize_quadtree(mesh15)["min_cell_width"])
    np.testing.assert_allclose(w15, w0, rtol=0.15, atol=1e-4)


@pytest.mark.parametrize("layout_key", ["2D-CN"])
def test_quadtree_backend_parity_with_grid_twisted(layout_key: str):
    twist_rad = float(np.deg2rad(15.0))
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    section = build_section_view(
        airfoil, resolve_system_type(layout_key), twist_angle_rad=twist_rad
    )
    airfoil_b = airfoil.rotate(twist_rad)
    grid = SDFGrid.from_airfoil(airfoil_b, nx=384, ny=180)
    cache = SectionEvalCache()
    grid_report = SectionPropertiesReport(
        section, grid, backend="grid", eval_cache=cache
    )
    qt_spec = quadtree_spec_for_airfoil(airfoil_b, chord=1.0, max_depth=8, max_cells=120000)
    qt_report = SectionPropertiesReport(
        section,
        grid,
        backend="quadtree",
        quadtree_spec=qt_spec,
        chord=1.0,
        twist_rad=twist_rad,
        airfoil_for_spec=airfoil_b,
    )
    for label in section.labels:
        g = grid_report[label]
        q = qt_report[label]
        if g["area"] < 1e-5:
            continue
        rtol = 0.12 if g["area"] < 0.002 else 0.08
        np.testing.assert_allclose(q["area"], g["area"], rtol=rtol, atol=5e-4)
        np.testing.assert_allclose(q["cx"], g["cx"], atol=0.05)
        np.testing.assert_allclose(q["cy"], g["cy"], atol=0.05)
