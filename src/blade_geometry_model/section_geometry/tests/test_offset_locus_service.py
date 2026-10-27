"""Offset-locus facade matches helpers and yields a mid-depth ring."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.section_geometry.constructive import (
    helpers as helpers_mod,
)
from blade_geometry_model.section_geometry.constructive import (
    offset_locus,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)


def test_airfoil_offset_ring_facade_matches_helpers():
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    depth = 0.003
    via_facade = offset_locus.airfoil_offset_ring(airfoil, depth)
    via_helpers = helpers_mod.airfoil_offset_ring(airfoil, depth)
    np.testing.assert_allclose(via_facade, via_helpers)


def test_mid_depth_locus_nonempty_naca2412():
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    skin = 0.006
    ring = offset_locus.locus_at_depth(airfoil, 0.5 * skin)
    assert ring.ndim == 2 and ring.shape[1] == 2
    assert ring.shape[0] >= 3
    assert np.isfinite(ring).all()
    cavity = offset_locus.cavity_polygon(airfoil, skin)
    assert not cavity.is_empty
    assert float(cavity.area) > 0.0