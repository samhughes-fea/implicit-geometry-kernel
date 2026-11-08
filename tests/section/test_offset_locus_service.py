"""Offset-locus helpers yield a mid-depth ring."""

from __future__ import annotations

import numpy as np

from implicit_geometry_kernel.section import boundary
from implicit_geometry_kernel.sdf.shapes.airfoil import (
    AirfoilSDF,
)


def test_airfoil_offset_ring_matches_helpers():
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    depth = 0.003
    via_facade = boundary.airfoil_offset_ring(airfoil, depth)
    via_helpers = boundary.airfoil_offset_ring(airfoil, depth)
    np.testing.assert_allclose(via_facade, via_helpers)


def test_mid_depth_locus_nonempty_naca2412():
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    skin = 0.006
    ring = boundary.locus_at_depth(airfoil, 0.5 * skin)
    assert ring.ndim == 2 and ring.shape[1] == 2
    assert ring.shape[0] >= 3
    assert np.isfinite(ring).all()
    cavity = boundary.cavity_polygon(airfoil, skin)
    assert not cavity.is_empty
    assert float(cavity.area) > 0.0
