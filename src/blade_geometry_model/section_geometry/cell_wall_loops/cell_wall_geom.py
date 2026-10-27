"""Skin-ring helpers for the void-cell wall loops: ordering, closure, area."""

from __future__ import annotations

from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from blade_geometry_model.section_contract.shell_midline_strip import ShellMidlineStrip


def shoelace_area(pts: NDArray[np.float64]) -> float:
    x = pts[:, 0]
    y = pts[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def drop_closing_vertex(pts: NDArray[np.float64], *, tol: float = 1e-12) -> NDArray[np.float64]:
    pts = np.asarray(pts, dtype=float)
    if pts.shape[0] > 2 and np.hypot(*(pts[0] - pts[-1])) < tol:
        return pts[:-1].copy()
    return pts.copy()


_SUCTION = "suction_surface"
_PRESSURE = "pressure_surface"
_LEADING = "leading_edge"
_TRAILING = "trailing_edge"
_BAY = "bay"

_SIDE_ALIASES = {
    "ss": _SUCTION,
    "upper": _SUCTION,
    _SUCTION: _SUCTION,
    "ps": _PRESSURE,
    "lower": _PRESSURE,
    _PRESSURE: _PRESSURE,
}
_ARC_ALIASES = {
    "le": _LEADING,
    _LEADING: _LEADING,
    "te": _TRAILING,
    _TRAILING: _TRAILING,
    "bay": _BAY,
}


def _normalise_skin_side(side: object | None) -> str | None:
    if side is None:
        return None
    key = str(side).strip().lower()
    return _SIDE_ALIASES.get(key)


def _normalise_skin_arc(arc: object | None) -> str | None:
    if arc is None:
        return None
    key = str(arc).strip().lower()
    return _ARC_ALIASES.get(key)


def _skin_ring_order_key(s: ShellMidlineStrip) -> tuple[int, int]:
    """Order hierarchical arcs: suction TE→LE, then pressure LE→TE."""
    side = _normalise_skin_side(s.skin_side)
    arc = _normalise_skin_arc(s.skin_arc)
    bay = int(s.bay_index) if s.bay_index is not None else 0
    side_rank = 0 if side == _SUCTION else 1 if side == _PRESSURE else 2
    if side == _SUCTION:
        # Half-ring runs trailing_edge → bay_i → leading_edge.
        if arc == _TRAILING:
            arc_rank = 0
        elif arc == _BAY:
            arc_rank = 1 + bay
        elif arc == _LEADING:
            arc_rank = 10_000
        else:
            arc_rank = 5_000
    else:
        # Half-ring runs leading_edge → bay_i → trailing_edge.
        if arc == _LEADING:
            arc_rank = 0
        elif arc == _BAY:
            arc_rank = 1 + bay
        elif arc == _TRAILING:
            arc_rank = 10_000
        else:
            arc_rank = 5_000
    return (side_rank, arc_rank)


def _append_oriented(
    poly: NDArray[np.float64],
    nxt: NDArray[np.float64],
    *,
    tol: float,
) -> NDArray[np.float64]:
    end = poly[-1]
    d0 = float(np.hypot(*(nxt[0] - end)))
    d1 = float(np.hypot(*(nxt[-1] - end)))
    if d1 < d0:
        nxt = nxt[::-1].copy()
        d0 = d1
    if d0 > max(tol, 1e-4):
        raise ValueError(f"Could not stitch skin strips into a closed ring (gap={d0:.3e} m)")
    nxt = drop_closing_vertex(nxt, tol=max(tol, 1e-9))
    if np.hypot(*(nxt[0] - poly[-1])) < max(tol, 1e-9):
        return np.vstack([poly, nxt[1:]])
    return np.vstack([poly, nxt])


def stitch_skin_ring(skin_strips: Sequence[ShellMidlineStrip], *, tol: float = 1e-9) -> NDArray[np.float64]:
    """Chain one or more skin strips into a single closed (open-vertex) ring.

    When strips carry hierarchical ``skin_side`` / ``skin_arc`` fields, prefer
    ring order suction (trailing→leading) then pressure (leading→trailing).
    Otherwise fall back to greedy tip chaining (legacy closed ring or
    unlabelled open segments). Short side/arc aliases are normalised at the
    edge of this helper.
    """
    if not skin_strips:
        raise ValueError("Expected at least one skin strip for cell wall loops")
    if len(skin_strips) == 1:
        return drop_closing_vertex(
            np.asarray(
                skin_strips[0].midline_twisted_pitch_axis_blade_frame, dtype=float
            ),
            tol=tol,
        )

    usable = [
        (i, s, np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float))
        for i, s in enumerate(skin_strips)
        if np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float).shape[0]
        >= 2
    ]
    if not usable:
        raise ValueError("Skin strips have no usable polylines")

    sides = [_normalise_skin_side(s.skin_side) for _i, s, _p in usable]
    hierarchical = all(side in (_SUCTION, _PRESSURE) for side in sides)
    if hierarchical and any(_normalise_skin_arc(s.skin_arc) for _i, s, _p in usable):
        ordered = sorted(usable, key=lambda item: _skin_ring_order_key(item[1]))
        try:
            poly = drop_closing_vertex(ordered[0][2].copy(), tol=tol)
            for _i, _s, pts in ordered[1:]:
                poly = _append_oriented(poly, pts.copy(), tol=tol)
            return drop_closing_vertex(poly, tol=max(tol, 1e-6))
        except ValueError:
            pass  # fall through to greedy tip chaining

    remaining = {i: p for i, _s, p in usable}

    def _tip_key(i: int) -> float:
        p = remaining[i]
        return float(max(p[0, 0], p[-1, 0]))

    start_i = max(remaining.keys(), key=_tip_key)
    poly = remaining.pop(start_i).copy()
    poly = drop_closing_vertex(poly, tol=tol)

    while remaining:
        end = poly[-1]
        best_i = None
        best_rev = False
        best_d = 1e300
        for i, p in remaining.items():
            d0 = float(np.hypot(*(p[0] - end)))
            d1 = float(np.hypot(*(p[-1] - end)))
            if d0 < best_d:
                best_d, best_i, best_rev = d0, i, False
            if d1 < best_d:
                best_d, best_i, best_rev = d1, i, True
        if best_i is None or best_d > max(tol, 1e-4):
            raise ValueError(
                f"Could not stitch skin strips into a closed ring (gap={best_d:.3e} m)"
            )
        nxt = remaining.pop(best_i)
        if best_rev:
            nxt = nxt[::-1].copy()
        nxt = drop_closing_vertex(nxt, tol=tol)
        if np.hypot(*(nxt[0] - poly[-1])) < max(tol, 1e-9):
            poly = np.vstack([poly, nxt[1:]])
        else:
            poly = np.vstack([poly, nxt])

    return drop_closing_vertex(poly, tol=max(tol, 1e-6))


__all__ = [
    "shoelace_area",
    "drop_closing_vertex",
    "stitch_skin_ring",
]
