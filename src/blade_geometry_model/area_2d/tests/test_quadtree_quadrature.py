"""Tests for quadtree section property quadrature."""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.props_emit.area_mesh import (
    _component_supports,
)
from blade_geometry_model.section_geometry.sdf.grid import SDFGrid
from blade_geometry_model.area_2d import QuadtreeSpec, build_shared_quadtree
from blade_geometry_model.area_2d.quadrature import (
    section_properties_from_quadtree,
)


def _circle_sampler(radius: float = 0.25):
    def sampler(pts: np.ndarray) -> np.ndarray:
        arr = np.asarray(pts, dtype=float).reshape(-1, 2)
        return np.hypot(arr[:, 0], arr[:, 1]) - radius

    return sampler


def test_quadtree_circle_area_matches_grid():
    spec = QuadtreeSpec(-0.6, 0.6, -0.6, 0.6, max_depth=7, min_cell_size=0.01, refine_band=0.02)
    sampler = _circle_sampler(0.25)
    mesh = build_shared_quadtree(spec, {"circle": sampler})
    qt_props = section_properties_from_quadtree(mesh, mesh.phi_by_label, "circle")

    grid = SDFGrid.from_bbox(-0.6, 0.6, -0.6, 0.6, nx=400, ny=400)
    phi = grid.eval(lambda x, y: np.hypot(x, y) - 0.25)
    grid_props = grid.section_properties(phi)

    np.testing.assert_allclose(qt_props["area"], grid_props["area"], rtol=0.02, atol=1e-3)
    np.testing.assert_allclose(qt_props["cx"], grid_props["cx"], atol=0.02)
    np.testing.assert_allclose(qt_props["cy"], grid_props["cy"], atol=0.02)


def test_quadtree_multi_field_refines_union():
    spec = QuadtreeSpec(-0.5, 0.5, -0.5, 0.5, max_depth=5, min_cell_size=0.03, refine_band=0.05)
    samplers = {"a": _circle_sampler(0.15), "b": _circle_sampler(0.35)}
    mesh = build_shared_quadtree(spec, samplers)
    leaves = mesh.leaf_cells()
    assert len(leaves) > 1
    single = build_shared_quadtree(spec, {"a": samplers["a"]})
    assert len(leaves) >= len(single.leaf_cells())


def test_shared_quadtree_leaf_cells_are_square():
    spec = QuadtreeSpec(-0.6, 0.6, -0.6, 0.6, max_depth=6, min_cell_size=0.02, refine_band=0.03)
    mesh = build_shared_quadtree(spec, {"circle": _circle_sampler(0.25)})
    leaves = mesh.leaf_cells()
    assert leaves
    for cell in leaves:
        np.testing.assert_allclose(cell.width, cell.height, rtol=0.0, atol=1e-12)


def test_rectangular_fallback_when_square_cells_disabled():
    spec = QuadtreeSpec(
        -0.6,
        0.6,
        -0.3,
        0.3,
        max_depth=4,
        min_cell_size=0.05,
        refine_band=0.05,
        square_cells=False,
    )
    mesh = build_shared_quadtree(spec, {"circle": _circle_sampler(0.2)})
    roots = mesh.roots
    assert len(roots) == 1
    assert roots[0].width > roots[0].height


def test_component_support_bbox_is_local():
    spec = QuadtreeSpec(-0.5, 0.5, -0.5, 0.5, min_cell_size=0.02)
    samplers = {"a": _circle_sampler(0.12)}
    supports = _component_supports(samplers, spec, n=32)
    xmin, xmax, ymin, ymax = supports["a"]
    assert xmax - xmin < 0.5
    assert ymax - ymin < 0.5


def test_effective_refine_band_scales_with_cell_size():
    from blade_geometry_model.area_2d.quadtree import _cell_needs_refine

    phi_near = np.array([0.01, 0.01, 0.01, 0.01, 0.01])
    assert _cell_needs_refine(phi_near, band_ceiling=0.02, band_scale=1.5, cell_diag=0.01)
    assert not _cell_needs_refine(phi_near, band_ceiling=0.02, band_scale=1.5, cell_diag=0.005)


def test_boundary_moments_batch_matches_per_cell():
    """Vectorized boundary moments must reproduce the per-cell marching-squares path."""
    from blade_geometry_model.area_2d.quadrature import (
        _boundary_moments_batch,
        _inside_polygon_from_cell,
        _polygon_moments,
    )
    from blade_geometry_model.area_2d.quadtree import QuadCell

    rng = np.random.default_rng(0)
    # Cover all 16 corner sign patterns plus random magnitudes, with assorted cells.
    bounds_rows = []
    phi_rows = []
    for code in range(16):
        signs = np.array([1.0 if (code >> i) & 1 else -1.0 for i in range(4)])
        for _ in range(40):
            x0 = rng.uniform(-1.0, 0.5)
            y0 = rng.uniform(-1.0, 0.5)
            x1 = x0 + rng.uniform(0.05, 0.6)
            y1 = y0 + rng.uniform(0.05, 0.6)
            # signs convention: inside = phi <= 0 ⇒ inside corner gets negative phi
            phi = -signs * rng.uniform(0.01, 1.0, size=4)
            bounds_rows.append([x0, x1, y0, y1])
            phi_rows.append(phi)
    bounds = np.array(bounds_rows, dtype=float)
    phi = np.array(phi_rows, dtype=float)

    # Reference: per-cell accumulation.
    ref = np.zeros(6)
    for b, pc in zip(bounds, phi):
        cell = QuadCell(b[0], b[1], b[2], b[3], level=0)
        poly = _inside_polygon_from_cell(cell, pc)
        if poly is None:
            continue
        ref += np.array(_polygon_moments(poly))

    got = _boundary_moments_batch(bounds, phi)
    np.testing.assert_allclose(got, ref, rtol=1e-12, atol=1e-12)
