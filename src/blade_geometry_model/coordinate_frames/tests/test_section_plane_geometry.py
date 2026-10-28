"""Unit tests for section-plane geometry and plot contract helpers."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from blade_geometry_model.coordinate_frames import (
    PIVOT_FRACTION_DEFAULT,
    naca_to_twisted_pitch_axis_blade_frame,
    naca_to_zero_twist_pitch_axis_station_frame,
    section_plot_axes,
    setup_section_axes,
    twisted_pitch_axis_blade_frame_to_naca,
)


def test_naca_to_station_leading_edge_at_negative_y():
    c = 1.0
    f = float(PIVOT_FRACTION_DEFAULT)
    verts = np.array([[0.0, 0.0], [f * c, 0.0], [c, 0.0]], dtype=float)
    pts = naca_to_zero_twist_pitch_axis_station_frame(verts, c)
    assert pts[0, 1] == pytest.approx(-f * c)
    assert pts[1, 1] == pytest.approx(0.0)
    assert pts[2, 1] == pytest.approx((1.0 - f) * c)
    assert pts[:, 0] == pytest.approx(0.0)


def test_naca_to_blade_yle_negative_yte_positive():
    c = 2.0
    f = float(PIVOT_FRACTION_DEFAULT)
    verts = np.array([[0.0, 0.1], [c, -0.05]], dtype=float)
    pts = naca_to_twisted_pitch_axis_blade_frame(verts, 0.0, c)
    assert pts[0, 1] < 0.0 < pts[1, 1]
    assert pts[0, 1] == pytest.approx(-f * c)
    back = twisted_pitch_axis_blade_frame_to_naca(pts, 0.0, c)
    assert back == pytest.approx(verts)


def test_section_plot_axes_order_and_invert_idempotent():
    pts = np.array([[0.1, -0.3], [0.0, 0.0], [-0.05, 0.5]], dtype=float)
    h, v = section_plot_axes(pts)
    assert h == pytest.approx(pts[:, 1])
    assert v == pytest.approx(pts[:, 0])

    fig, ax = plt.subplots()
    try:
        ax.plot(h, v)
        setup_section_axes(ax)
        assert ax.xaxis_inverted()
        setup_section_axes(ax)
        assert ax.xaxis_inverted()
    finally:
        plt.close(fig)
