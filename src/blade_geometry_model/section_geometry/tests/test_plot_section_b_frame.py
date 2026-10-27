"""plot_section must place the pitch pivot at the B-frame origin."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from blade_geometry_model.coordinate_frames import (
    PIVOT_FRACTION_DEFAULT,
    naca_to_twisted_pitch_axis_blade_frame,
)
from blade_geometry_model.section_geometry.sections.multicell import (
    MultiCellSection,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.sdf.grid import SDFGrid
from blade_geometry_model.section_geometry.visualisation.plot_section import (
    plot_section,
    rotated_naca_grid_to_section_plot_axes,
)


@pytest.fixture(scope="module")
def af():
    return AirfoilSDF.from_naca("2412", chord=1.0)


def test_rotated_naca_grid_zero_twist_shifts_pitch_pivot_to_origin(af):
    c = float(af.chord)
    f = float(PIVOT_FRACTION_DEFAULT)
    X = np.array([[0.0, f * c, c]], dtype=float)
    Y = np.array([[0.0, 0.0, 0.0]], dtype=float)
    y_b, x_b = rotated_naca_grid_to_section_plot_axes(X, Y, twist_rad=0.0, chord_m=c)
    # LE / pitch / TE on the chordline → (Y_B, X_B) = (-f c, 0), (0, 0), ((1-f)c, 0)
    assert y_b[0, 0] == pytest.approx(-f * c)
    assert y_b[0, 1] == pytest.approx(0.0)
    assert y_b[0, 2] == pytest.approx((1.0 - f) * c)
    assert x_b == pytest.approx(0.0)


def test_rotated_naca_grid_matches_naca_to_blade(af):
    c = float(af.chord)
    twist = float(np.deg2rad(12.0))
    verts = np.asarray(af.vertices, dtype=float)

    y0, x0 = rotated_naca_grid_to_section_plot_axes(
        verts[:, 0], verts[:, 1], twist_rad=0.0, chord_m=c
    )
    pts0 = naca_to_twisted_pitch_axis_blade_frame(verts, 0.0, c)
    assert y0 == pytest.approx(pts0[:, 1])
    assert x0 == pytest.approx(pts0[:, 0])

    # Synthetic LE-rotated cloud of the same verts (MultiCellSection twist space).
    ct, st = float(np.cos(twist)), float(np.sin(twist))
    qx = ct * verts[:, 0] - st * verts[:, 1]
    qy = st * verts[:, 0] + ct * verts[:, 1]
    y_r, x_r = rotated_naca_grid_to_section_plot_axes(qx, qy, twist_rad=twist, chord_m=c)
    pts_b = naca_to_twisted_pitch_axis_blade_frame(verts, twist, c)
    assert y_r == pytest.approx(pts_b[:, 1], abs=1e-12)
    assert x_r == pytest.approx(pts_b[:, 0], abs=1e-12)


def test_plot_section_and_sdf_and_mesh_xlim_straddle_pitch(af):
    from blade_geometry_model.section_geometry.visualisation.plot_sdf import (
        plot_sdf_field,
    )
    from blade_geometry_model.section_geometry.visualisation.plot_section import (
        plot_section_from_mesh,
    )
    from blade_geometry_model.section_geometry.props_emit.export import (
        SectionPropertiesReport,
    )
    from blade_geometry_model.section_geometry.props_emit.props_backend import (
        quadtree_spec_for_airfoil,
    )

    twin = MultiCellSection.twin_web(af, twist_angle=0.0)
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=80, ny=40)
    f = float(PIVOT_FRACTION_DEFAULT)
    c = float(af.chord)

    fig, ax = plot_section(twin, grid, show_legend=False, title="t")
    try:
        xmin, xmax = ax.get_xlim()
        assert min(xmin, xmax) < 0.0 < max(xmin, xmax)
        assert min(xmin, xmax) == pytest.approx(-f * c, abs=0.15)
    finally:
        plt.close(fig)

    phi = grid.eval(twin.airfoil)
    fig, ax = plot_sdf_field(phi, grid, chord_m=c, twist_rad=0.0, title="phi")
    try:
        xmin, xmax = ax.get_xlim()
        assert min(xmin, xmax) < 0.0 < max(xmin, xmax)
        assert ax.xaxis_inverted()
    finally:
        plt.close(fig)

    spec = quadtree_spec_for_airfoil(af, chord=c, max_depth=5, max_cells=50_000)
    rep = SectionPropertiesReport(
        twin, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=af
    )
    fig, ax = plot_section_from_mesh(
        rep.mesh, show_legend=False, title="m", chord_m=c, section=twin
    )
    try:
        xmin, xmax = ax.get_xlim()
        assert min(xmin, xmax) < 0.0 < max(xmin, xmax)
        assert ax.xaxis_inverted()
    finally:
        plt.close(fig)


def test_plot_section_xlim_straddles_pitch_origin(af):
    twin = MultiCellSection.twin_web(af, twist_angle=0.0)
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=120, ny=60)
    fig, ax = plot_section(twin, grid, show_legend=False, title="t")
    try:
        xmin, xmax = ax.get_xlim()
        # Inverted axis (LE on RHS): data still has Y_LE < 0 < Y_TE.
        assert min(xmin, xmax) < 0.0 < max(xmin, xmax)
        # Pitch fraction ~0.333c from LE → LE near -0.333, TE near +0.667.
        f = float(PIVOT_FRACTION_DEFAULT)
        assert min(xmin, xmax) == pytest.approx(-f * float(af.chord), abs=0.12)
        assert max(xmin, xmax) == pytest.approx((1.0 - f) * float(af.chord), abs=0.12)
    finally:
        plt.close(fig)


def test_plot_section_twisted_skin_bbox_matches_blade_verts(af):
    """Nonzero twist must B-resample; remapped NACA grids clip the section."""
    twist = float(np.deg2rad(15.0))
    twin = MultiCellSection.twin_web(af, twist_angle=twist)
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=100, ny=50)
    c = float(af.chord)
    verts_b = naca_to_twisted_pitch_axis_blade_frame(
        np.asarray(af.vertices, dtype=float), twist, c
    )
    y_lo, y_hi = float(verts_b[:, 1].min()), float(verts_b[:, 1].max())
    x_lo, x_hi = float(verts_b[:, 0].min()), float(verts_b[:, 0].max())

    fig, ax = plot_section(twin, grid, show_legend=False, title="t15", auto_limits=True)
    try:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        assert min(xmin, xmax) < 0.0 < max(xmin, xmax)
        assert ax.xaxis_inverted()
        # Occupancy window should cover most of the twisted airfoil B bbox
        # (not a thin clipped band from a warped NACA lattice).
        assert min(xmin, xmax) <= y_lo + 0.08 * c
        assert max(xmin, xmax) >= y_hi - 0.08 * c
        assert ymin <= x_lo + 0.05 * c
        assert ymax >= x_hi - 0.05 * c
        # Vertical span must be meaningful (warped remap collapses flap extent).
        assert (ymax - ymin) >= 0.55 * (x_hi - x_lo)
    finally:
        plt.close(fig)


def test_blade_plot_sampling_grid_rectilinear_and_query_roundtrip(af):
    from blade_geometry_model.section_geometry.visualisation.plot_section import (
        blade_plot_sampling_grid,
        section_query_from_blade,
    )

    twist = float(np.deg2rad(15.0))
    c = float(af.chord)
    Yh, Xv, QX, QY = blade_plot_sampling_grid(
        af, twist_rad=twist, chord_m=c, nx=40, ny=20
    )
    assert Yh.shape == Xv.shape == QX.shape == QY.shape
    # meshgrid(yb, xb): Yh[i,j]=yb[j] (const along column), Xv[i,j]=xb[i] (const along row).
    assert np.allclose(Yh[:, 5], Yh[0, 5])
    assert np.allclose(Xv[3, :], Xv[3, 0])
    assert float(Yh.min()) < 0.0 < float(Yh.max())
    pts_b = np.column_stack([Xv.ravel(), Yh.ravel()])
    qx2, qy2 = section_query_from_blade(pts_b, twist_rad=twist, chord_m=c)
    assert qx2 == pytest.approx(QX.ravel(), abs=1e-12)
    assert qy2 == pytest.approx(QY.ravel(), abs=1e-12)
