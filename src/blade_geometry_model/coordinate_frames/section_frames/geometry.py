"""NACA ↔ zero-twist pitch-axis station ↔ twisted pitch-axis blade transforms.

Short aliases ``naca_to_S``, ``rotate_S_to_B``, ``chord_to_blade``, … remain for
callers; prefer the long primary names in new code.
"""

from __future__ import annotations

import numpy as np

PIVOT_FRACTION_DEFAULT = 0.333
"""Pitch-axis chordline as a fraction of chord, measured from the LE.

Matches ``blade_geometry.pivot_fraction``. Prefer passing the configured value
through; this default exists so the transform is usable standalone.
"""


def naca_to_zero_twist_pitch_axis_station_frame(
    verts_naca: np.ndarray,
    chord_m: float,
    pivot_fraction: float = PIVOT_FRACTION_DEFAULT,
) -> np.ndarray:
    """Raw NACA ``(x_chord, y_thickness)`` → canonical S ``(X_S, Y_S)``.

    .. math::
        X_S = y_{\\mathrm{thickness}}, \\qquad
        Y_S = x_{\\mathrm{chord}} - f_{\\mathrm{pivot}}\\, c

    Both steps matter. The **swap** puts chord-normal first, matching
    \\(X\\) = flap; the **pivot shift** moves the origin from the LE onto the
    pitch-axis chordline, which is what makes \\(Y_{LE} < 0 < Y_{TE}\\).

    Parameters
    ----------
    verts_naca : np.ndarray
        ``(N, 2)`` raw NACA vertices, LE at ``x = 0``.
    chord_m : float
        Section chord [m].
    pivot_fraction : float
        Pitch-axis position as a fraction of chord from the LE.
    """
    v = np.asarray(verts_naca, dtype=np.float64)
    if v.ndim != 2 or v.shape[1] != 2:
        raise ValueError(f"verts_naca must be (N, 2); got {v.shape}")
    return np.column_stack([v[:, 1], v[:, 0] - float(pivot_fraction) * float(chord_m)])


def zero_twist_pitch_axis_station_frame_to_naca(
    pts_S: np.ndarray,
    chord_m: float,
    pivot_fraction: float = PIVOT_FRACTION_DEFAULT,
) -> np.ndarray:
    """Inverse of :func:`naca_to_zero_twist_pitch_axis_station_frame`."""
    p = np.asarray(pts_S, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 2:
        raise ValueError(f"pts_S must be (N, 2); got {p.shape}")
    return np.column_stack([p[:, 1] + float(pivot_fraction) * float(chord_m), p[:, 0]])


def rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame(
    pts_S: np.ndarray, twist_rad: float
) -> np.ndarray:
    """Rotate canonical S points by twist → B (twisted pitch-axis blade frame).

    .. math::
        X_B = \\cos\\theta\\, X_S + \\sin\\theta\\, Y_S, \\qquad
        Y_B = -\\sin\\theta\\, X_S + \\cos\\theta\\, Y_S

    S is the zero-twist section view; B carries the geometric twist.

    .. warning::
       **This is ``R(-theta)``, and that is convention, not a typo.** On a
       canonical ``(X, Y)`` array the S → B rotation is ``R(-theta)`` — as used
       by ``visualisation/demos/plot_structural_frames.py`` and
       ``feasibility/visualisation/plots.py``, the two places already working in
       canonical order. A *swapped* ``(Y, X)`` array needs ``R(+theta)`` for the
       same physical rotation, because ``P R(+t) P = R(-t)`` for the swap ``P``.
       Passing a swapped array here rotates it the wrong way.
    """
    p = np.asarray(pts_S, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 2:
        raise ValueError(f"pts_S must be (N, 2); got {p.shape}")
    t = float(twist_rad)
    if abs(t) < 1e-15:
        return p.copy()
    c, s = np.cos(t), np.sin(t)
    return np.column_stack([c * p[:, 0] + s * p[:, 1], -s * p[:, 0] + c * p[:, 1]])


def rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame_swapped_order(
    pts_S: np.ndarray, twist_rad: float
) -> np.ndarray:
    """S → B for a legacy **swapped** ``(Y, X)`` section-plane array.

    Exactly ``P · rotate_...(P · pts)`` for the component swap ``P``.
    """
    p = np.asarray(pts_S, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 2:
        raise ValueError(f"pts_S must be (N, 2); got {p.shape}")
    return rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame(
        p[:, ::-1], twist_rad
    )[:, ::-1]


def rotate_twisted_pitch_axis_blade_frame_to_zero_twist_pitch_axis_station_frame(
    pts_B: np.ndarray, twist_rad: float
) -> np.ndarray:
    """Inverse of S→B rotation (rotation by ``-twist_rad``)."""
    return rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame(
        pts_B, -float(twist_rad)
    )


def naca_to_twisted_pitch_axis_blade_frame(
    verts_naca: np.ndarray,
    twist_rad: float,
    chord_m: float,
    pivot_fraction: float = PIVOT_FRACTION_DEFAULT,
) -> np.ndarray:
    """Raw NACA ``(x_chord, y_thickness)`` → canonical B ``(X_B, Y_B)``.

    Full handoff: NACA → S then S → B. Use at the geometry→shell boundary.
    """
    return rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame(
        naca_to_zero_twist_pitch_axis_station_frame(
            verts_naca, chord_m, pivot_fraction
        ),
        twist_rad,
    )


def twisted_pitch_axis_blade_frame_to_naca(
    pts_B: np.ndarray,
    twist_rad: float,
    chord_m: float,
    pivot_fraction: float = PIVOT_FRACTION_DEFAULT,
) -> np.ndarray:
    """Inverse of :func:`naca_to_twisted_pitch_axis_blade_frame`."""
    return zero_twist_pitch_axis_station_frame_to_naca(
        rotate_twisted_pitch_axis_blade_frame_to_zero_twist_pitch_axis_station_frame(
            pts_B, twist_rad
        ),
        chord_m,
        pivot_fraction,
    )


# --- Short aliases (S / B letter names) ---------------------------------------
naca_to_S = naca_to_zero_twist_pitch_axis_station_frame
S_to_naca = zero_twist_pitch_axis_station_frame_to_naca
rotate_S_to_B = (
    rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame
)
rotate_S_to_B_swapped_order = (
    rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame_swapped_order
)
rotate_B_to_S = (
    rotate_twisted_pitch_axis_blade_frame_to_zero_twist_pitch_axis_station_frame
)
chord_to_blade = naca_to_twisted_pitch_axis_blade_frame
blade_to_chord = twisted_pitch_axis_blade_frame_to_naca
