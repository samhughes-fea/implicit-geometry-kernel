"""Regression: product shell channel map + LE-on-RHS S mapping."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.coordinate_frames import (
    blade_resultants_to_shell_channels,
    naca_to_zero_twist_pitch_axis_station_frame,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF


def test_blade_resultants_to_shell_channels_flap_on_Mz():
    ch = blade_resultants_to_shell_channels(
        M_flapwise=1000.0,
        M_edgewise=200.0,
        V_flapwise=300.0,
        V_edgewise=40.0,
    )
    assert ch.Mz == 1000.0
    assert ch.My == 200.0
    assert ch.Vz == 300.0
    assert ch.Vy == 40.0


def test_naca_to_station_leading_edge_negative_y():
    """RH section: +Y toward TE so LE sits at negative Y (RHS after plot invert)."""
    af = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    poly_S = naca_to_zero_twist_pitch_axis_station_frame(af.vertices, chord_m=1.0)
    y_le = -0.333  # -pivot * chord with LE at x_naca=0
    y_te = float(np.max(poly_S[:, 1]))
    assert y_le < y_te
    assert float(np.min(poly_S[:, 1])) <= y_le + 1e-9
