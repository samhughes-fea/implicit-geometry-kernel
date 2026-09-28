"""Snap continuous spar-cap heights onto a ply pitch (BL-019).

Geometry stays material-agnostic: callers pass `pitch_m` (catalogue
`PlyGrid.pitch` for the spar-cap stack). After cavity feasibility reduces a
cap height off-grid, flooring onto the pitch makes midline `thickness_m`
match the ply-count ceil used at Catalog FI bind.
"""

from __future__ import annotations

import math

_ONE_PLY_REL_TOL = 1e-9


def snap_cap_height_at_most(height_m: float, pitch_m: float | None) -> float:
    """Greatest on-grid height that does not exceed `height_m`.

    Returns `height_m` unchanged when `pitch_m` is absent/invalid.
    When not even one ply fits (`height_m < pitch`), returns ``0.0`` — geometric
    run-out, matching ``PlyGrid.largest_at_most``. Leaving a sub-ply continuous
    strip made the ply-count ceil (1 ply) disagree with geometry and broke Catalog FI.
    """
    if pitch_m is None:
        return float(height_m)
    p = float(pitch_m)
    c = float(height_m)
    if not (p > 0.0) or not math.isfinite(p) or not math.isfinite(c):
        return float(height_m)
    if c < p:
        return 0.0
    return math.floor(c / p + _ONE_PLY_REL_TOL) * p


def snap_cap_heights(
    h_upper: float,
    h_lower: float,
    pitch_m: float | None,
) -> tuple[float, float]:
    """Snap a suction/pressure pair after feasibility clamp."""
    return (
        snap_cap_height_at_most(h_upper, pitch_m),
        snap_cap_height_at_most(h_lower, pitch_m),
    )


__all__ = ["snap_cap_height_at_most", "snap_cap_heights"]
