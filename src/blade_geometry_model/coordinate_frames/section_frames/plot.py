"""Section-plane plot contract (B/S/S′ axes, labels, +x3 normal)."""

from __future__ import annotations

import numpy as np


SECTION_PLOT_XLABEL = r"$Y_B$ (edgewise) [m]"
SECTION_PLOT_YLABEL = r"$X_B$ (flapwise) [m]"


def section_plot_axes(pts_B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Canonical B section points → ``(horizontal, vertical)`` plot arrays.

    Section plots draw **horizontal =** \\(Y_B\\) (edgewise) and **vertical =**
    \\(X_B\\) (flapwise), i.e. the transpose of the stored ``(X_B, Y_B)`` order.
    Pair with :func:`setup_section_axes`, which applies the labels and the
    leading-edge-on-RHS inversion.

    Plotting ``pts[:, 0]`` horizontally renders the section **rotated 90°**.
    """
    p = np.asarray(pts_B, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 2:
        raise ValueError(f"pts_B must be (N, 2); got {p.shape}")
    return p[:, 1], p[:, 0]


def setup_section_axes(
    ax,
    *,
    xlabel: str | None = None,
    ylabel: str | None = None,
    invert_for_le_on_rhs: bool = True,
    equal_aspect: bool = True,
) -> None:
    """Apply the canonical section-plot axis convention to ``ax``.

    Labels the axes \\(Y_B\\) / \\(X_B\\) and inverts the horizontal axis so the
    leading edge appears on the RHS.

    .. important::
       The inversion is **display-only**. Stored coordinates keep
       \\(+Y_B\\) toward the TE (``Y_LE < 0 < Y_TE``) — never emulate the flip by
       negating data, which silently breaks every downstream geometric test.

    For aerofoil section figures, call :func:`fit_square_section_axes` after
    data are plotted so the data window and axes box are square.

    Idempotent per axis: repeated calls will not double-invert.
    """
    ax.set_xlabel(SECTION_PLOT_XLABEL if xlabel is None else xlabel)
    ax.set_ylabel(SECTION_PLOT_YLABEL if ylabel is None else ylabel)
    if equal_aspect:
        ax.set_aspect("equal", adjustable="datalim")
    if invert_for_le_on_rhs and not ax.xaxis_inverted():
        ax.invert_xaxis()


def fit_square_section_axes(
    ax,
    *,
    pad_frac: float = 0.08,
    box_aspect: float | None = 1.0,
) -> None:
    """Expand current data limits to a square window (for aerofoil section views).

    Preserves horizontal inversion. When ``box_aspect`` is set (default ``1``),
    the axes *box* in the figure is also square so rotated sections are not
    squeezed into a wide chord-length panel.

    Call after data are plotted and after :func:`setup_section_axes`.
    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    x0, x1 = (min(xmin, xmax), max(xmin, xmax))
    y0, y1 = (min(ymin, ymax), max(ymin, ymax))
    cx = 0.5 * (x0 + x1)
    cy = 0.5 * (y0 + y1)
    half = 0.5 * max(x1 - x0, y1 - y0, 1e-9)
    half *= 1.0 + max(float(pad_frac), 0.0)
    if xmin > xmax:  # already inverted (LE on RHS)
        ax.set_xlim(cx + half, cx - half)
    else:
        ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    if box_aspect is not None:
        ax.set_box_aspect(float(box_aspect))


def plot_section_polyline(ax, pts_B: np.ndarray, /, **kwargs):
    """``ax.plot`` a canonical-B polyline with the section axis convention applied.

    Convenience over :func:`section_plot_axes` for the common case; the caller
    still owns styling via ``kwargs``.
    """
    h, v = section_plot_axes(pts_B)
    return ax.plot(h, v, **kwargs)


def section_normal_2d(t: np.ndarray) -> np.ndarray:
    """Unit normal 90° CCW of a section-plane tangent ``t``.

    Canonical owner of the section-plane ``+x3`` convention; the shell mesh's
    ``panel_x3_orientation.normal_2d`` delegates here. With a canonical
    ``(X_B, Y_B)`` array this is ``Z_B × t``, making ``(x1, x2, x3)``
    right-handed for ``x1 = +Z_B``.
    """
    v = np.asarray(t, dtype=np.float64).ravel()[:2]
    n = float(np.linalg.norm(v))
    if n < 1e-14:
        return np.array([1.0, 0.0], dtype=np.float64)
    v = v / n
    return np.array([-float(v[1]), float(v[0])], dtype=np.float64)

