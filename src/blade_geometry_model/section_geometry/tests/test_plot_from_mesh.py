"""Quadtree-mesh plot renderer: area tie-out and structural checks.

`plot_section_from_mesh` renders section fills directly from the shared quadtree
mesh (reusing the property backend's `mesh.phi_by_label`) instead of a second
dense SDF grid. The fill is the union of interior-leaf rectangles + boundary-leaf
marching-squares clips, so it must tile the same φ<0 region the property quadrature
integrates. These tests tie the rendered fill area to the integrated section area
and confirm the renderer emits fill + boundary collections per component.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.collections import LineCollection, PolyCollection

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import quadtree_spec_for_airfoil
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.props_emit.export import SectionPropertiesReport
from blade_geometry_model.section_geometry.visualisation.plot import plot_section_from_mesh
from blade_geometry_model.area_2d.quadrature import (
    label_boundary_segments,
    label_fill_polygons,
)


def _polygon_area(p: np.ndarray) -> float:
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))


def _report(key: str, twist_deg: float = 0.0):
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    layout = resolve_system_type(key)
    section = build_section_view(
        airfoil, layout, twist_angle_rad=float(np.deg2rad(twist_deg))
    )
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=6, max_cells=150000)
    return SectionPropertiesReport(
        section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=airfoil
    )


@pytest.mark.parametrize("key,twist", [("2D-CN", 0.0), ("3C-F", 12.0), ("0B", 0.0)])
def test_fill_polygon_area_matches_property_area(key, twist):
    rep = _report(key, twist)
    mesh = rep.mesh
    assert mesh is not None
    for label in mesh.phi_by_label:
        fill_area = sum(_polygon_area(p) for p in label_fill_polygons(mesh, label))
        prop_area = float(rep[label]["area"])
        # Fill is the same interior-rect + boundary-clip partition the quadrature
        # integrates, so areas agree to round-off. SS half-annulus camber clips can
        # leave a tiny hollow-island residual (~5e-6 m²) vs property quadrature.
        assert fill_area == pytest.approx(prop_area, rel=2e-5, abs=5e-6), label


def test_renderer_emits_fill_and_boundary_collections():
    rep = _report("2D-CN")
    mesh = rep.mesh
    fig, ax = plot_section_from_mesh(mesh, title="t", chord_m=1.0)
    try:
        polys = [c for c in ax.collections if isinstance(c, PolyCollection)]
        lines = [c for c in ax.collections if isinstance(c, LineCollection)]
        n_labels = len(mesh.phi_by_label)
        # One fill + one boundary collection per non-empty label.
        assert 0 < len(polys) <= n_labels
        assert 0 < len(lines) <= n_labels
        # Canonical B: pitch pivot at origin, Y_LE < 0 < Y_TE, LE on RHS after invert.
        xmin, xmax = ax.get_xlim()
        assert min(xmin, xmax) < 0.0 < max(xmin, xmax)
        assert ax.xaxis_inverted()
        # Mesh domain was NACA LE-origin; remapped xlim must not equal raw mesh.spec.
        assert ax.get_xlim() != pytest.approx((mesh.spec.x_max, mesh.spec.x_min))
    finally:
        plt.close(fig)


def test_boundary_segments_are_pairs():
    rep = _report("2D-CN")
    mesh = rep.mesh
    for label in mesh.phi_by_label:
        segs = label_boundary_segments(mesh, label)
        assert segs.ndim == 3 and segs.shape[1:] == (2, 2)
