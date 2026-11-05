"""Contour helpers for 2-D area mesh zero-level extraction."""

from __future__ import annotations

import numpy as np


def ensure_closed(points: np.ndarray) -> np.ndarray:
    """Ensure the polygon is closed (last point == first point)."""
    if len(points) == 0:
        return points
    if np.allclose(points[0], points[-1]):
        return points
    return np.vstack([points, points[0]])


def smooth_contour(points: np.ndarray, passes: int = 3) -> np.ndarray:
    """1-4-1 Laplacian smoothing on a closed contour."""
    pts = ensure_closed(points)
    core = pts[:-1].copy()
    n = len(core)
    if n < 5:
        return pts
    for _ in range(passes):
        core = (
            0.25 * np.roll(core, 1, axis=0)
            + 0.5 * core
            + 0.25 * np.roll(core, -1, axis=0)
        )
    return ensure_closed(core)
