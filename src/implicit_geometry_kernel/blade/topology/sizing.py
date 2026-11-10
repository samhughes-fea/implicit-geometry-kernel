"""In-plane size field ``h`` on skin-mid / web polylines (CSG viz first).

Pinned keys stay put. Counts follow ``h``; this module does not author ``n_chord``.
"""
from __future__ import annotations

import numpy as np

H_MAX_OVER_C = 0.04
EPS_OVER_C = 0.0015
H_MIN_OVER_T = 0.125
KAPPA_FLOOR = 1e-12


def station_live_thickness(cut) -> float:
    """Thinnest live skin or web gauge at the station (caps are partitions)."""
    ts = [
        float(w.thickness)
        for w in cut.walls
        if w.kind in ("skin", "web") and float(w.thickness) > 0.0
    ]
    return min(ts) if ts else 0.0


def size_field_h(kappa, chord: float, t: float) -> np.ndarray:
    """``h = max(h_min, min(h_max, sagitta))`` at each sample.

    ``t/4`` is the wall-meta quadtree ``min_cell``, not a midline spacing cap.
    Putting it in ``min(...)`` forces ~1 mm nodes and swamps sagitta / ``h_max``.
    """
    kappa = np.asarray(kappa, float)
    c = max(float(chord), 1e-12)
    t = max(float(t), 0.0)
    h_max = H_MAX_OVER_C * c
    h_sagitta = np.sqrt(8.0 * (EPS_OVER_C * c) / np.maximum(kappa, KAPPA_FLOOR))
    h_min = t * H_MIN_OVER_T if t > 0.0 else 0.0
    return np.maximum(h_min, np.minimum(h_max, h_sagitta))


def polyline_kappa(pts: np.ndarray) -> np.ndarray:
    """Discrete curvature ``|dθ|/ds`` at vertices of an open polyline."""
    xy = np.asarray(pts, float)[:, :2]
    n = len(xy)
    kap = np.zeros(n)
    if n < 3:
        return kap
    dxy = np.diff(xy, axis=0)
    seg = np.linalg.norm(dxy, axis=1)
    tang = dxy / np.maximum(seg[:, None], 1e-15)
    cross = tang[:-1, 0] * tang[1:, 1] - tang[:-1, 1] * tang[1:, 0]
    dot = np.clip(np.sum(tang[:-1] * tang[1:], axis=1), -1.0, 1.0)
    dtheta = np.arctan2(cross, dot)
    ds = 0.5 * (seg[:-1] + seg[1:])
    kap[1:-1] = np.abs(dtheta) / np.maximum(ds, 1e-15)
    kap[0] = kap[1]
    kap[-1] = kap[-2]
    return kap


def resample_open_h(pts: np.ndarray, chord: float, t: float) -> np.ndarray:
    """Arc-length resample with local spacing ``≈ h``. Endpoints are kept."""
    pts = np.asarray(pts, float)
    if len(pts) == 0:
        return np.zeros((2, 3))
    if pts.ndim != 2:
        raise ValueError("points must be (n, d)")
    if pts.shape[1] == 2:
        pts = np.column_stack([pts, np.zeros(len(pts))])
    if len(pts) == 1:
        return np.repeat(pts[:1], 2, axis=0)
    seg = np.linalg.norm(np.diff(pts[:, :2], axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 1e-15:
        out = np.repeat(pts[:1], 2, axis=0)
        out[-1] = pts[-1]
        return out
    h = size_field_h(polyline_kappa(pts), chord, t)
    inv = 1.0 / np.maximum(h, 1e-15)
    weight = np.concatenate([[0.0], np.cumsum(0.5 * (inv[:-1] + inv[1:]) * seg)])
    n = max(2, int(np.ceil(float(weight[-1]) - 1e-12)) + 1)
    u = np.linspace(0.0, float(weight[-1]), n)
    out = np.empty((n, pts.shape[1]))
    for d in range(pts.shape[1]):
        out[:, d] = np.interp(u, weight, pts[:, d])
    out[0] = pts[0]
    out[-1] = pts[-1]
    return out
