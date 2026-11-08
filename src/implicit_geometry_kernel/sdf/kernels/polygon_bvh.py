"""Cached edges and AABB for lofted polygon SDF (same phi as the edge loop)."""

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


def build_polygon_accel(verts: np.ndarray, nx: int = 16, ny: int = 16) -> dict:
    """AABB plus a uniform grid of edge indices for batched lofted queries."""
    verts = np.ascontiguousarray(verts[:, :2], dtype=np.float64)
    n = int(verts.shape[0])
    ax = np.ascontiguousarray(verts[:, 0], dtype=np.float64)
    ay = np.ascontiguousarray(verts[:, 1], dtype=np.float64)
    bx = np.empty(n, dtype=np.float64)
    by = np.empty(n, dtype=np.float64)
    bx[:-1] = ax[1:]
    by[:-1] = ay[1:]
    bx[-1] = ax[0]
    by[-1] = ay[0]
    xmin = float(min(ax.min(), bx.min()))
    xmax = float(max(ax.max(), bx.max()))
    ymin = float(min(ay.min(), by.min()))
    ymax = float(max(ay.max(), by.max()))
    pad = 1e-12 * max(1.0, xmax - xmin, ymax - ymin)
    dx = (xmax - xmin + 2.0 * pad) / float(nx) if nx else 1.0
    dy = (ymax - ymin + 2.0 * pad) / float(ny) if ny else 1.0
    buckets: list[list[int]] = [[] for _ in range(int(nx) * int(ny))]
    x0, y0 = xmin - pad, ymin - pad
    for e in range(n):
        ex0, ex1 = (ax[e], bx[e]) if ax[e] <= bx[e] else (bx[e], ax[e])
        ey0, ey1 = (ay[e], by[e]) if ay[e] <= by[e] else (by[e], ay[e])
        gx0 = max(0, min(nx - 1, int((ex0 - x0) / dx)))
        gx1 = max(0, min(nx - 1, int((ex1 - x0) / dx)))
        gy0 = max(0, min(ny - 1, int((ey0 - y0) / dy)))
        gy1 = max(0, min(ny - 1, int((ey1 - y0) / dy)))
        for gy in range(gy0, gy1 + 1):
            for gx in range(gx0, gx1 + 1):
                buckets[gy * nx + gx].append(e)
    counts = np.array([len(b) for b in buckets], dtype=np.int32)
    offsets = np.zeros(int(nx) * int(ny) + 1, dtype=np.int32)
    offsets[1:] = np.cumsum(counts)
    flat = np.empty(int(offsets[-1]), dtype=np.int32)
    for i, bucket in enumerate(buckets):
        if bucket:
            flat[offsets[i] : offsets[i + 1]] = np.asarray(bucket, dtype=np.int32)
    return {
        "ax": ax,
        "ay": ay,
        "bx": bx,
        "by": by,
        "nx": int(nx),
        "ny": int(ny),
        "offsets": offsets,
        "flat": flat,
        "aabb": (xmin, xmax, ymin, ymax),
    }
