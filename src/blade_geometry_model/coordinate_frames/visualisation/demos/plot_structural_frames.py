"""Plot B / S / S′ structural frames on a NACA section (LE on RHS).

Canon: ``coordinate_frames/section_frames/README.md`` — storage
``(X_*, Y_*)`` with ``+Y`` toward TE (``Y_LE < 0 < Y_TE``); plot horizontal
``Y_*``, vertical ``X_*``, invert horizontal (display-only) so LE is on the RHS.

Run::

    python -m blade_geometry_model.coordinate_frames.visualisation.demos.plot_structural_frames
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from blade_geometry_model.coordinate_frames import (
    PIVOT_FRACTION_DEFAULT,
    fit_square_section_axes,
    naca_to_zero_twist_pitch_axis_station_frame,
    rotate_S_to_B,
    setup_section_axes,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF

# This demo deliberately imports the *production* transforms rather than
# defining its own. A demo that reimplements the maths can agree with the docs
# while the product disagrees with both -- which is exactly what happened here:
# naca_to_zero_twist_pitch_axis_station_frame (then aliased as naca_to_S) lived
# only in this file, and production applied a bare rotation.
_OUT = Path(__file__).resolve().parent.parent / "output"
_PIVOT_FRAC = PIVOT_FRACTION_DEFAULT
_CHORD_M = 1.0
_TWIST_DEG = 15.0
_RR = 0.50  # optional R-frame context


def S_to_B(points_S: NDArray[np.float64], twist_deg: float) -> NDArray[np.float64]:
    """Degrees-taking wrapper over the canonical ``rotate_S_to_B``.

    Kept only so the plotting code below can pass degrees; the transform itself
    lives in ``coordinate_frames``.
    """
    return rotate_S_to_B(points_S, np.deg2rad(float(twist_deg)))


def _centroid_xy(poly: NDArray[np.float64]) -> tuple[float, float]:
    x, y = poly[:, 0], poly[:, 1]
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    cross = x * y2 - x2 * y
    a = 0.5 * float(np.sum(cross))
    if abs(a) < 1e-18:
        return float(np.mean(x)), float(np.mean(y))
    cx = float(np.sum((x + x2) * cross) / (6.0 * a))
    cy = float(np.sum((y + y2) * cross) / (6.0 * a))
    return cx, cy


def _draw_triad(ax, origin, ex, ey, labels, scale: float) -> None:
    o = np.asarray(origin, dtype=float)
    for vec, lab, color in (
        (np.asarray(ex, dtype=float), labels[0], "C3"),
        (np.asarray(ey, dtype=float), labels[1], "C0"),
    ):
        v = vec / (np.linalg.norm(vec) + 1e-30) * scale
        ax.annotate(
            "",
            xy=o + v,
            xytext=o,
            arrowprops=dict(arrowstyle="->", color=color, lw=2),
        )
        ax.text(*(o + 1.12 * v), lab, color=color, fontsize=11, fontweight="bold")


def _panel(ax, poly, title: str, xlab: str, ylab: str, triad, le_xy, origin, r_note: str | None) -> None:
    # Plot coords: horizontal = Y, vertical = X (flap up / edge horizontal).
    ax.fill(poly[:, 1], poly[:, 0], facecolor="0.85", edgecolor="k", lw=1.0)
    ax.plot(poly[:, 1], poly[:, 0], "k-", lw=1.0)
    ax.plot([le_xy[1]], [le_xy[0]], "ro", ms=8, label="LE")
    ax.plot([origin[1]], [origin[0]], "ks", ms=6, label="origin")
    scale = 0.22 * float(np.ptp(poly[:, 1]))
    # Triad in plot space: +X = up, +Y = increasing Y (displays left after invert).
    _draw_triad(ax, (origin[1], origin[0]), (0.0, 1.0), (1.0, 0.0), triad, scale)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    if r_note:
        ax.text(0.98, 0.02, r_note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
    # LE on RHS + equal aspect, then square data window / axes box.
    setup_section_axes(ax, xlabel=xlab, ylabel=ylab)
    fit_square_section_axes(ax)


def main() -> int:
    _OUT.mkdir(parents=True, exist_ok=True)
    af = AirfoilSDF.from_naca("2412", n_points=240, chord=_CHORD_M)
    poly_S = naca_to_zero_twist_pitch_axis_station_frame(af.vertices, _CHORD_M)
    poly_B = S_to_B(poly_S, _TWIST_DEG)
    cx_s, cy_s = _centroid_xy(poly_S)
    poly_Sp = np.column_stack([poly_S[:, 0] - cx_s, poly_S[:, 1] - cy_s])

    # NACA LE at x_chord=0 → Y_S = -pivot * chord (LHS in data; RHS after invert)
    le_S = np.array([0.0, -_PIVOT_FRAC * _CHORD_M], dtype=np.float64)
    le_B = S_to_B(le_S.reshape(1, 2), _TWIST_DEG)[0]
    le_Sp = np.array([le_S[0] - cx_s, le_S[1] - cy_s])

    r_note = f"R context: r/R = {_RR:.2f}"
    # Wide enough for three square axes boxes (box_aspect=1); avoid bbox_inches='tight'.
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.2))

    _panel(
        axes[0],
        poly_B,
        f"B frame (twist = {_TWIST_DEG:.0f}°)\n$X_B$ flap, $Y_B$ edge, $Z_B$ span",
        r"$Y_B$ (edgewise) → +Y toward TE; LE on RHS",
        r"$X_B$ (flapwise)",
        (r"$+X_B$", r"$+Y_B$"),
        le_B,
        (0.0, 0.0),
        r_note,
    )
    _panel(
        axes[1],
        poly_S,
        f"S frame (zero-twist)\n$X_S$ chord-normal, $Y_S$ chord-parallel, $Z_S$ span",
        r"$Y_S$ (chord-parallel) → +Y toward TE; LE on RHS",
        r"$X_S$ (chord-normal)",
        (r"$+X_S$", r"$+Y_S$"),
        le_S,
        (0.0, 0.0),
        r_note,
    )
    _panel(
        axes[2],
        poly_Sp,
        "S′ frame (centroid origin)\nsame axes as S",
        r"$Y_{S'}$ → +Y toward TE; LE on RHS",
        r"$X_{S'}$",
        (r"$+X_{S'}$", r"$+Y_{S'}$"),
        le_Sp,
        (0.0, 0.0),
        "pivot (S) marked relative shift in docs",
    )
    # mark S pivot in S′ panel (plot coords = (Y, X)); re-square after the extra artist
    axes[2].plot([-cy_s], [-cx_s], "b^", ms=7, label="S pivot")
    axes[2].legend(loc="upper left", fontsize=8)
    fit_square_section_axes(axes[2])

    fig.suptitle("Structural frames (canonical) — NACA 2412", fontsize=13)
    fig.tight_layout()
    out = _OUT / "structural_frames_B_S_Sprime.png"
    # Do not use bbox_inches='tight' — it fights square box_aspect panels.
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    # LE at negative Y (data); appears on RHS after invert_xaxis
    y_te = float(np.max(poly_S[:, 1]))
    y_le = float(le_S[1])
    if not (y_le < y_te):
        raise RuntimeError(f"LE-at-negative-Y check failed: Y_LE={y_le}, Y_TE={y_te}")
    print(f"LE-at-negative-Y ok: Y_LE={y_le:.4f} < Y_TE={y_te:.4f} (LE on RHS after invert)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
