"""Cached edge arrays and box of a lofted station polygon (same phi as the edge loop)."""

from __future__ import annotations

import numpy as np

from implicit_geometry_kernel.sdf.kernels.polygon_numba import _sdf_polygon_numba_core


def sdf_polygon_bvh(xf: np.ndarray, yf: np.ndarray, accel: dict) -> np.ndarray:
    """Evaluate using cached edges; distances match ``sdf_polygon``."""
    return _sdf_polygon_numba_core(
        np.ascontiguousarray(xf, dtype=np.float64),
        np.ascontiguousarray(yf, dtype=np.float64),
        accel["ax"],
        accel["ay"],
        accel["bx"],
        accel["by"],
    )


def build_polygon_accel(verts: np.ndarray) -> dict:
    """The station polygon's edge arrays and axis-aligned box, cached per station.

    ``sdf_polygon_bvh`` runs the compiled edge loop over these arrays, so ``phi`` is
    bit-for-bit that of ``sdf_polygon``; the saving is that the arrays are built once
    per station rather than once per query. No spatial search structure is built:
    an earlier version filled a uniform grid of edge buckets that no query read
    (SD-004).
    """
    verts = np.ascontiguousarray(verts[:, :2], dtype=np.float64)
    ax = np.ascontiguousarray(verts[:, 0], dtype=np.float64)
    ay = np.ascontiguousarray(verts[:, 1], dtype=np.float64)
    bx = np.ascontiguousarray(np.roll(ax, -1))
    by = np.ascontiguousarray(np.roll(ay, -1))
    aabb = (float(np.min(ax)), float(np.max(ax)), float(np.min(ay)), float(np.max(ay)))
    return {"ax": ax, "ay": ay, "bx": bx, "by": by, "aabb": aabb}
