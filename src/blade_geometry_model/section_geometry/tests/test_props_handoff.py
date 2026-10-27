"""Geometry → section_properties peer handoffs (fiber extents / topology)."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    fiber_extents_from_airfoil,
    fiber_extents_from_outline,
    outer_polygon_from_airfoil,
    topology_for_section_properties,
)


def test_fiber_extents_from_outline_about_centroid() -> None:
    outline = np.array([[0.0, -2.0], [4.0, -2.0], [4.0, 2.0], [0.0, 2.0]], dtype=float)
    ext = fiber_extents_from_outline(outline, Cx=2.0, Cy=0.0)
    assert ext is not None
    assert np.isclose(ext["x_max"], 2.0)
    assert np.isclose(ext["y_max"], 2.0)


def test_fiber_extents_from_airfoil_and_topology() -> None:
    airfoil = AirfoilSDF.from_naca("0012", n_points=80, chord=1.0)
    ext = fiber_extents_from_airfoil(airfoil, Cx=0.5, Cy=0.0)
    assert ext is not None
    assert ext["x_max"] > 0.0 and ext["y_max"] > 0.0
    poly = outer_polygon_from_airfoil(airfoil)
    assert poly is not None and poly.area > 0.0
    topo = topology_for_section_properties(airfoil=airfoil)
    assert "outer_polygon" in topo
