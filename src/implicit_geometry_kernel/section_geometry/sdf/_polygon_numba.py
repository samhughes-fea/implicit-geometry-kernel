"""Numba edge-loop SDF for closed polygons (required production path).

Imported from ``primitives.sdf_polygon``. Numba is a hard dependency of the
package; there is no optional / silent fallback for production evaluation.
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _sdf_polygon_numba_core(
    xf: np.ndarray,
    yf: np.ndarray,
    ax: np.ndarray,
    ay: np.ndarray,
    bx: np.ndarray,
    by: np.ndarray,
) -> np.ndarray:
    m = xf.shape[0]
    n = ax.shape[0]
    out = np.empty(m, dtype=np.float64)
    for i in range(m):
        px = xf[i]
        py = yf[i]
        d_min = 1.0e300
        wn = 0
        for e in range(n):
            eax = ax[e]
            eay = ay[e]
            ebx = bx[e]
            eby = by[e]
            abx = ebx - eax
            aby = eby - eay
            apx = px - eax
            apy = py - eay
            denom = abx * abx + aby * aby + 1e-30
            t = (apx * abx + apy * aby) / denom
            if t < 0.0:
                t = 0.0
            elif t > 1.0:
                t = 1.0
            dx = apx - t * abx
            dy = apy - t * aby
            d2 = dx * dx + dy * dy
            if d2 < d_min:
                d_min = d2
            c3 = abx * (py - eay) - aby * (px - eax)
            if eay <= py:
                if eby > py and c3 > 0.0:
                    wn += 1
            else:
                if eby <= py and c3 < 0.0:
                    wn -= 1
        sign = 1.0 if wn == 0 else -1.0
        out[i] = sign * np.sqrt(d_min)
    return out


def sdf_polygon_numba(xf: np.ndarray, yf: np.ndarray, verts: np.ndarray) -> np.ndarray:
    """Evaluate polygon SDF with a Numba edge loop (1-D point arrays)."""
    n = int(verts.shape[0])
    ax = np.ascontiguousarray(verts[:, 0], dtype=np.float64)
    ay = np.ascontiguousarray(verts[:, 1], dtype=np.float64)
    bx = np.empty(n, dtype=np.float64)
    by = np.empty(n, dtype=np.float64)
    bx[:-1] = ax[1:]
    by[:-1] = ay[1:]
    bx[-1] = ax[0]
    by[-1] = ay[0]
    return _sdf_polygon_numba_core(
        np.ascontiguousarray(xf, dtype=np.float64),
        np.ascontiguousarray(yf, dtype=np.float64),
        ax,
        ay,
        bx,
        by,
    )
