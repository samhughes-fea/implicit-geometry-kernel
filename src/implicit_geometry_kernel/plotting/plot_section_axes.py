"""Shared axes, colour, and NACA-to-B helpers for section plots."""

from __future__ import annotations

import numpy as np

from .palette import REGION_FALLBACK, region_color

_DEFAULT_COLOR = REGION_FALLBACK
_DEFAULT_CORE_ALPHA = 0.22


def fit_square_section_axes(
    ax,
    *,
    pad_frac: float = 0.08,
    box_aspect: float | None = 1.0,
) -> None:
    """Expand current data limits to a square window (for twisted B-frame views).

    Canon owner: ``coordinate_frames.section_frames.plot.fit_square_section_axes``.
    This is a compatibility delegate.
    """
    from .frames import fit_square_section_axes as _fit_square_section_axes

    _fit_square_section_axes(ax, pad_frac=pad_frac, box_aspect=box_aspect)


def apply_structural_section_axes(
    ax,
    *,
    invert_horizontal: bool = True,
    xlabel: str | None = None,
    ylabel: str | None = None,
) -> None:
    """Finish a section plot in the canonical structural frame.

    Axes are **horizontal =** ``Y_B`` (edgewise), **vertical =** ``X_B``
    (flapwise), with the horizontal axis inverted so the LE (negative ``Y_B``)
    appears on the RHS. The inversion is display-only.

    .. important::
       **This sets up the axes; it does not transpose your data.** Which
       component to plot horizontally depends on the frame the data is in:

       * **Canonical B** (``midline_twisted_pitch_axis_blade_frame``,
         ``polyline_twisted_pitch_axis_blade_frame``,
         ``leading_edge_twisted_pitch_axis_blade_frame``) is
         stored ``(X_B, Y_B)`` — plot ``pts[:, 1]`` horizontally. Use
         ``coordinate_frames.section_plot_axes`` rather than indexing by hand.
       * **SDF grids / rotated-NACA solids** are authored with LE at
         ``x_chord = 0``. Convert to B with
         :func:`rotated_naca_grid_to_section_plot_axes` (pivot shift + twist)
         before labeling axes as ``Y_B`` / ``X_B``.

       Do **not** overlay raw NACA polylines on a B-labeled plot without
       ``coordinate_frames.naca_to_twisted_pitch_axis_blade_frame``.

    See ``geometry_utilities/coordinate_frames/README.md`` (source of truth).

    Idempotent w.r.t. horizontal inversion: safe to call after ``set_xlim``.
    """
    from .frames import setup_section_axes

    setup_section_axes(
        ax,
        xlabel=xlabel,
        ylabel=ylabel,
        invert_for_le_on_rhs=bool(invert_horizontal),
    )
    if not invert_horizontal and ax.xaxis_inverted():
        ax.invert_xaxis()


def rotated_naca_grid_to_section_plot_axes(
    X_naca: np.ndarray,
    Y_naca: np.ndarray,
    *,
    twist_rad: float = 0.0,
    chord_m: float,
    pivot_fraction: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Map SDF-grid ``(x_chord, y_thickness)`` → plot ``(Y_B, X_B)``.

    Preferred for **zero twist** (pure pitch-pivot shift). For nonzero twist,
    prefer :func:`blade_plot_sampling_grid` + SDF query in B — remapping a
    chord-aligned NACA mesh into B warps the sample lattice and clips fills.

    At zero twist: ``Y_B = x_chord - f_pivot * c``, ``X_B = y_thickness``.
    """
    from implicit_geometry_kernel.frames import (
        PIVOT_FRACTION_DEFAULT,
        naca_to_twisted_pitch_axis_blade_frame,
    )

    X = np.asarray(X_naca, dtype=float)
    Y = np.asarray(Y_naca, dtype=float)
    if X.shape != Y.shape:
        raise ValueError(f"X_naca and Y_naca shapes must match; got {X.shape} vs {Y.shape}")
    c_m = float(chord_m)
    if not np.isfinite(c_m) or c_m <= 0.0:
        raise ValueError(f"chord_m must be positive and finite; got {chord_m!r}")
    f_pivot = (
        PIVOT_FRACTION_DEFAULT if pivot_fraction is None else float(pivot_fraction)
    )
    t = float(twist_rad)
    qx = X.ravel()
    qy = Y.ravel()
    if abs(t) < 1e-15:
        naca = np.column_stack([qx, qy])
    else:
        # Inverse of MultiCellSection's rotate_field about the LE: R(-twist).
        ct, st = float(np.cos(t)), float(np.sin(t))
        naca = np.column_stack([ct * qx + st * qy, -st * qx + ct * qy])
    pts_b = naca_to_twisted_pitch_axis_blade_frame(naca, t, c_m, f_pivot)
    y_b = pts_b[:, 1].reshape(X.shape)
    x_b = pts_b[:, 0].reshape(X.shape)
    return y_b, x_b


def section_query_from_blade(
    pts_b: np.ndarray,
    *,
    twist_rad: float,
    chord_m: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Map canonical B ``(X_B, Y_B)`` → MultiCellSection SDF query coords.

    Twisted solids live in rotated-NACA (LE-pivot) space, not full
    ``naca_to_twisted_pitch_axis_blade_frame`` B. Convert B → NACA via
    ``twisted_pitch_axis_blade_frame_to_naca``, then apply ``R(+twist)`` about
    the LE when twisted.
    """
    from implicit_geometry_kernel.frames import (
        twisted_pitch_axis_blade_frame_to_naca,
    )

    naca = twisted_pitch_axis_blade_frame_to_naca(
        pts_b, float(twist_rad), float(chord_m)
    )
    t = float(twist_rad)
    if abs(t) < 1e-12:
        return naca[:, 0], naca[:, 1]
    ct, st = float(np.cos(t)), float(np.sin(t))
    qx = ct * naca[:, 0] - st * naca[:, 1]
    qy = st * naca[:, 0] + ct * naca[:, 1]
    return qx, qy


def blade_plot_sampling_grid(
    airfoil,
    *,
    twist_rad: float,
    chord_m: float | None = None,
    nx: int = 400,
    ny: int = 200,
    pad_frac: float = 0.08,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build a rectilinear B-frame sample mesh for SDF contour / pcolor plots.

    Returns
    -------
    Y_plot, X_plot, QX, QY
        ``Y_plot``/``X_plot`` are meshgrids for ``contourf(Y_plot, X_plot, …)``
        (horizontal = ``Y_B``, vertical = ``X_B``). ``QX``/``QY`` are the matching
        rotated-NACA query coordinates for ``section[label](QX, QY)``.
    """
    from implicit_geometry_kernel.frames import (
        naca_to_twisted_pitch_axis_blade_frame,
    )

    c_m = float(chord_m) if chord_m is not None else float(
        getattr(airfoil, "chord", 0.0) or 0.0
    )
    if not np.isfinite(c_m) or c_m <= 0.0:
        raise ValueError("blade_plot_sampling_grid requires a positive chord")
    verts = np.asarray(airfoil.vertices, dtype=float)
    verts_b = naca_to_twisted_pitch_axis_blade_frame(verts, float(twist_rad), c_m)
    span = float(max(np.ptp(verts_b[:, 0]), np.ptp(verts_b[:, 1]), 1e-9))
    pad = float(pad_frac) * span
    xb = np.linspace(float(verts_b[:, 0].min()) - pad, float(verts_b[:, 0].max()) + pad, int(ny))
    yb = np.linspace(float(verts_b[:, 1].min()) - pad, float(verts_b[:, 1].max()) + pad, int(nx))
    # indexing='xy': Y_plot[i,j] varies with j (horizontal), X_plot with i (vertical)
    Y_plot, X_plot = np.meshgrid(yb, xb)
    pts_b = np.column_stack([X_plot.ravel(), Y_plot.ravel()])
    qx, qy = section_query_from_blade(
        pts_b, twist_rad=float(twist_rad), chord_m=c_m
    )
    return Y_plot, X_plot, qx.reshape(X_plot.shape), qy.reshape(X_plot.shape)


def _resolve_section_twist_chord(
    section=None,
    mesh=None,
    *,
    twist_rad: float | None = None,
    chord_m: float | None = None,
) -> tuple[float, float]:
    """Resolve twist [rad] and chord [m] for NACA→B plot remaps."""
    twist = float(twist_rad) if twist_rad is not None else float(
        getattr(section, "twist_rad", 0.0) or 0.0
    )
    c = float(chord_m) if chord_m is not None else 0.0
    if not (np.isfinite(c) and c > 0.0) and section is not None:
        af = getattr(section, "airfoil", None) or getattr(section, "_af", None)
        c = float(getattr(af, "chord", 0.0) or 0.0)
    if not (np.isfinite(c) and c > 0.0) and mesh is not None:
        spec = getattr(mesh, "spec", None)
        if spec is not None:
            # Quadtree domain is NACA-like; chord ≈ LE→TE span on x.
            c = float(getattr(spec, "x_max", 0.0) - getattr(spec, "x_min", 0.0))
    if not (np.isfinite(c) and c > 0.0):
        raise ValueError(
            "chord_m is required for B-frame section plots "
            "(pass chord_m=... or attach airfoil.chord)"
        )
    return twist, c


def _naca_xy_to_plot_hv(
    xy: np.ndarray,
    *,
    twist_rad: float,
    chord_m: float,
) -> np.ndarray:
    """``(N, 2)`` rotated-NACA → plot coords ``(N, 2)`` as ``(Y_B, X_B)``."""
    p = np.asarray(xy, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2:
        raise ValueError(f"xy must be (N, 2); got {p.shape}")
    h, v = rotated_naca_grid_to_section_plot_axes(
        p[:, 0], p[:, 1], twist_rad=float(twist_rad), chord_m=float(chord_m)
    )
    return np.column_stack([np.asarray(h, dtype=float), np.asarray(v, dtype=float)])


def _component_color(label):
    """Resolve a component label to the shared deck palette."""
    return region_color(str(getattr(label, "value", label)))


