"""Regression: LE-on-RHS S mapping."""

from __future__ import annotations

import numpy as np

from implicit_geometry_kernel.frames import naca_to_zero_twist_pitch_axis_station_frame
from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF


def test_naca_to_station_leading_edge_negative_y():
    """RH section: +Y toward TE so LE sits at negative Y (RHS after plot invert)."""
    af = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    poly_S = naca_to_zero_twist_pitch_axis_station_frame(af.vertices, chord_m=1.0)
    y_le = -0.333  # -pivot * chord with LE at x_naca=0
    y_te = float(np.max(poly_S[:, 1]))
    assert y_le < y_te
    assert float(np.min(poly_S[:, 1])) <= y_le + 1e-9
