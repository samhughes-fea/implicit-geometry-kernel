"""
Exact signed distance functions for 2-D geometric primitives.

All functions accept scalar or NumPy array inputs for (x, y) and return
an array of the same broadcast shape.

Sign convention:
    phi < 0  →  interior
    phi = 0  →  boundary
    phi > 0  →  exterior
"""

import time

import numpy as np

from implicit_geometry_kernel.profiling import profiling_enabled, record_sdf_polygon
from implicit_geometry_kernel.sdf.kernels.polygon_numba import sdf_polygon_numba


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _as_arrays(*args):
    return tuple(_ensure_f64(a) for a in args)


def _ensure_f64(a):
    """Return float64 ndarray with a fast-path for already-correct arrays."""
    if isinstance(a, np.ndarray) and a.dtype == np.float64:
        return a
    return np.asarray(a, dtype=float)


def _clamp(val, lo, hi):
    return np.clip(val, lo, hi)


def _length(vx, vy):
    return np.sqrt(vx**2 + vy**2)


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------

def sdf_circle(x, y, cx=0.0, cy=0.0, r=1.0):
    """Exact SDF for a circle."""
    x, y = _as_arrays(x, y)
    return _length(x - cx, y - cy) - r


def sdf_box(x, y, cx=0.0, cy=0.0, half_w=1.0, half_h=0.5):
    """Exact SDF for an axis-aligned rectangle."""
    x, y = _as_arrays(x, y)
    dx = np.abs(x - cx) - half_w
    dy = np.abs(y - cy) - half_h
    outside = _length(np.maximum(dx, 0.0), np.maximum(dy, 0.0))
    inside  = np.minimum(np.maximum(dx, dy), 0.0)
    return outside + inside


def sdf_half_plane(x, y, nx=0.0, ny=1.0, d=0.0):
    """Signed distance to a half-plane  n·p + d ≤ 0.

    The normal (nx, ny) points *outward* (into the positive / exterior region).
    """
    x, y = _as_arrays(x, y)
    n_mag = np.sqrt(nx**2 + ny**2)
    if n_mag < 1e-14:
        raise ValueError("Normal vector must be non-zero.")
    return (nx * x + ny * y + d) / n_mag


def sdf_segment(x, y, ax, ay, bx, by):
    """Exact SDF to a line segment AB (unsigned distance)."""
    x, y = _as_arrays(x, y)
    abx, aby = bx - ax, by - ay
    apx, apy = x - ax, y - ay
    t = _clamp((apx * abx + apy * aby) / (abx**2 + aby**2 + 1e-30), 0.0, 1.0)
    return _length(apx - t * abx, apy - t * aby)


def sdf_capsule(x, y, ax, ay, bx, by, r):
    """Exact SDF for a capsule (Minkowski sum of segment and disk of radius r)."""
    return sdf_segment(x, y, ax, ay, bx, by) - r


def sdf_rounded_box(x, y, cx=0.0, cy=0.0, half_w=1.0, half_h=0.5, r=0.1):
    """Exact SDF for a rectangle with uniformly rounded corners."""
    r = min(r, min(half_w, half_h))
    x, y = _as_arrays(x, y)
    dx = np.abs(x - cx) - half_w + r
    dy = np.abs(y - cy) - half_h + r
    outside = _length(np.maximum(dx, 0.0), np.maximum(dy, 0.0))
    inside  = np.minimum(np.maximum(dx, dy), 0.0)
    return outside + inside - r


def sdf_ellipse(x, y, cx=0.0, cy=0.0, rx=1.0, ry=0.5):
    """Exact signed distance to an axis-aligned ellipse (R4.10).

    The nearest point of the ellipse to ``(p0, p1)`` in the first quadrant, with
    semi-axes ``e0 >= e1``, is ``x_i = e_i^2 p_i / (s + e_i^2)``, where ``s`` is the
    root of ``F(s) = sum (e_i p_i / (s + e_i^2))^2 - 1``. ``F`` decreases on
    ``s > -e1^2`` and its root lies in ``[-e1^2 + e1 p1, -e1^2 + |(e0 p0, e1 p1)|]``,
    so bisection converges to it (Eberly, 2013). Points on an axis are solved in
    closed form.
    """
    x, y = _as_arrays(x, y)
    a = np.abs(np.asarray(x, dtype=float) - cx)
    b = np.abs(np.asarray(y, dtype=float) - cy)
    rx, ry = float(rx), float(ry)
    if rx >= ry:
        e0, e1, p0, p1 = rx, ry, a, b
    else:
        e0, e1, p0, p1 = ry, rx, b, a
    p0, p1 = np.broadcast_arrays(p0, p1)
    p0 = np.array(p0, dtype=float, copy=True)
    p1 = np.array(p1, dtype=float, copy=True)
    dist = np.empty_like(p0)

    general = (p0 > 0.0) & (p1 > 0.0)
    if np.any(general):
        q0, q1 = p0[general], p1[general]
        z0, z1 = e0 * q0, e1 * q1
        lo = -e1 * e1 + z1
        hi = -e1 * e1 + np.hypot(z0, z1)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            f = (z0 / (mid + e0 * e0)) ** 2 + (z1 / (mid + e1 * e1)) ** 2 - 1.0
            pos = f > 0.0
            lo = np.where(pos, mid, lo)
            hi = np.where(pos, hi, mid)
            if np.all(hi - lo <= 4.0 * np.finfo(float).eps * np.maximum(1.0, np.abs(mid))):
                break
        s = 0.5 * (lo + hi)
        n0 = e0 * e0 * q0 / (s + e0 * e0)
        n1 = e1 * e1 * q1 / (s + e1 * e1)
        dist[general] = np.hypot(q0 - n0, q1 - n1)

    on_minor = (p0 <= 0.0) & ~general
    dist[on_minor] = np.abs(p1[on_minor] - e1)

    on_major = (p1 <= 0.0) & (p0 > 0.0)
    if np.any(on_major):
        q0 = p0[on_major]
        if e0 > e1:
            reach = (e0 * e0 - e1 * e1) / e0
            inner = q0 < reach
            n0 = np.where(inner, e0 * e0 * q0 / (e0 * e0 - e1 * e1), e0)
            n1 = np.where(inner, e1 * np.sqrt(np.clip(1.0 - (n0 / e0) ** 2, 0.0, None)), 0.0)
            dist[on_major] = np.hypot(q0 - n0, n1)
        else:
            dist[on_major] = np.abs(q0 - e0)

    inside = (a / rx) ** 2 + (b / ry) ** 2 <= 1.0
    return np.where(inside, -dist, dist)


def sdf_oriented_box(x, y, cx=0.0, cy=0.0, half_w=1.0, half_h=0.5, angle=0.0):
    """Exact SDF for a rotated axis-aligned rectangle."""
    x, y = _as_arrays(x, y)
    cos_a = np.cos(-angle)
    sin_a = np.sin(-angle)
    xr = x - cx
    yr = y - cy
    xp = cos_a * xr - sin_a * yr
    yp = sin_a * xr + cos_a * yr
    dx = np.abs(xp) - half_w
    dy = np.abs(yp) - half_h
    outside = _length(np.maximum(dx, 0.0), np.maximum(dy, 0.0))
    inside  = np.minimum(np.maximum(dx, dy), 0.0)
    return outside + inside


def sdf_polygon(x, y, vertices):
    """Exact SDF for a closed polygon (winding-number sign).

    Uses the required Numba edge-loop kernel. There is no NumPy fallback on
    the production path; install ``numba`` (core package dependency).

    Parameters
    ----------
    vertices : array-like, shape (N, 2)
        Polygon vertices in order (need not be closed; last edge wraps).
    """
    x_b, y_b = np.broadcast_arrays(_ensure_f64(x), _ensure_f64(y))
    verts = _ensure_f64(vertices)
    if int(len(verts)) < 3:
        raise ValueError("Polygon must contain at least 3 vertices.")
    m = int(x_b.size)
    t0 = time.perf_counter()
    out = sdf_polygon_numba(x_b.ravel(), y_b.ravel(), verts).reshape(x_b.shape)
    elapsed = time.perf_counter() - t0
    if profiling_enabled():
        record_sdf_polygon(m, elapsed_s=elapsed)
    return out
