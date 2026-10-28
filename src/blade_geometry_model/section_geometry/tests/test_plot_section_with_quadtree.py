"""Tests for plot_section_with_quadtree."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.patches import Polygon, Rectangle

from blade_geometry_model.section_geometry.api.factory import build_section_view, resolve_system_type
from blade_geometry_model.section_geometry.props_emit.props_backend import quadtree_spec_for_airfoil
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.props_emit.export import SectionPropertiesReport
from blade_geometry_model.section_geometry.visualisation.plot import plot_section_with_quadtree


def _report(key: str = "2D-CN"):
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    layout = resolve_system_type(key)
    section = build_section_view(airfoil, layout, twist_angle_rad=0.0)
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=6, max_cells=150000)
    return section, SectionPropertiesReport(
        section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=airfoil
    )


def test_plot_section_with_quadtree_two_axes_and_rectangles():
    section, rep = _report()
    mesh = rep.mesh
    assert mesh is not None
    fig, axes = plot_section_with_quadtree(mesh, title="t", section=section, chord_m=1.0)
    try:
        assert len(fig.axes) == 2
        cells = [
            p for p in axes[1].patches
            if isinstance(p, (Rectangle, Polygon))
        ]
        assert len(cells) > 0
        xmin, xmax = axes[1].get_xlim()
        assert min(xmin, xmax) < 0.0 < max(xmin, xmax)
        assert axes[1].xaxis_inverted()
    finally:
        plt.close(fig)
