"""
Exact signed distance functions for 2-D geometric primitives.

All functions accept scalar or NumPy array inputs for (x, y) and return
an array of the same broadcast shape.

Sign convention:
    phi < 0  →  interior
    phi = 0  →  boundary
    phi > 0  →  exterior
"""

import numpy as np


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
    """Approximate SDF for an axis-aligned ellipse (Quilez iterative, 4 Newton steps)."""
    x, y = _as_arrays(x, y)
    px = np.abs(x - cx)
    py = np.abs(y - cy)
    tx = np.full_like(px, 0.707)
    ty = np.full_like(py, 0.707)
    for _ in range(4):
        ex = rx * tx
        ey = ry * ty
        qx = px - ex
        qy = py - ey
        r_ = _length(qx, qy)
        r_ = np.where(r_ < 1e-14, 1e-14, r_)
        lx = (ex - px) / (rx * rx)
        ly = (ey - py) / (ry * ry)
        l_ = _length(lx, ly)
        l_ = np.where(l_ < 1e-14, 1e-14, l_)
        tx = _clamp(lx / l_, 0.0, 1.0)
        ty = _clamp(ly / l_, 0.0, 1.0)
        t_mag = _length(tx, ty)
        t_mag = np.where(t_mag < 1e-14, 1e-14, t_mag)
        tx /= t_mag
        ty /= t_mag
    nearest_x = rx * tx
    nearest_y = ry * ty
    dist = _length(px - nearest_x, py - nearest_y)
    inside_mask = (px / rx)**2 + (py / ry)**2 <= 1.0
    return np.where(inside_mask, -dist, dist)


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


def sdf_polygon_numpy(x, y, vertices):
    """Exact SDF for a closed polygon (NumPy **test reference** only).

    Production code must use :func:`sdf_polygon` (Numba). Keep this for
    parity tests that compare against the vectorized reference.
    """
    x, y = np.broadcast_arrays(_ensure_f64(x), _ensure_f64(y))
    verts = _ensure_f64(vertices)
    n = int(len(verts))
    if n < 3:
        raise ValueError("Polygon must contain at least 3 vertices.")

    xf = x.ravel()
    yf = y.ravel()
    m = int(xf.size)

    ax = verts[:, 0]
    ay = verts[:, 1]
    bx = verts[(np.arange(n) + 1) % n, 0]
    by = verts[(np.arange(n) + 1) % n, 1]
    abx = bx - ax
    aby = by - ay
    ab_denom = (abx**2 + aby**2 + 1e-30)[:, None]

    d_out = np.empty(m, dtype=np.float64)
    winding_out = np.empty(m, dtype=np.int32)

    chunk_size = 1 << 16
    for start in range(0, m, chunk_size):
        stop = min(start + chunk_size, m)
        px = xf[start:stop][None, :]
        py = yf[start:stop][None, :]

        apx = px - ax[:, None]
        apy = py - ay[:, None]
        t = _clamp((apx * abx[:, None] + apy * aby[:, None]) / ab_denom, 0.0, 1.0)
        dx = apx - t * abx[:, None]
        dy = apy - t * aby[:, None]
        d_out[start:stop] = np.min(dx**2 + dy**2, axis=0)

        c1 = py >= ay[:, None]
        c2 = py < by[:, None]
        c3 = abx[:, None] * (py - ay[:, None]) - aby[:, None] * (px - ax[:, None])
        wn_up = np.sum(c1 & c2 & (c3 > 0.0), axis=0, dtype=np.int32)
        wn_dn = np.sum((~c1) & (~c2) & (c3 < 0.0), axis=0, dtype=np.int32)
        winding_out[start:stop] = wn_up - wn_dn

    sign = np.where(winding_out == 0, 1.0, -1.0)
    return (sign * np.sqrt(d_out)).reshape(x.shape)


def sdf_polygon(x, y, vertices):
    """Exact SDF for a closed polygon (winding-number sign).

    Uses the required Numba edge-loop kernel. There is no NumPy fallback on
    the production path; install ``numba`` (core package dependency).

    Parameters
    ----------
    vertices : array-like, shape (N, 2)
        Polygon vertices in order (need not be closed; last edge wraps).
    """
    import time

    from blade_geometry_model.section_geometry.sdf._polygon_numba import (
        sdf_polygon_numba,
    )

    x_b, y_b = np.broadcast_arrays(_ensure_f64(x), _ensure_f64(y))
    verts = _ensure_f64(vertices)
    if int(len(verts)) < 3:
        raise ValueError("Polygon must contain at least 3 vertices.")
    m = int(x_b.size)
    t0 = time.perf_counter()
    out = sdf_polygon_numba(x_b.ravel(), y_b.ravel(), verts).reshape(x_b.shape)
    elapsed = time.perf_counter() - t0
    try:
        from blade_geometry_model.section_geometry.sdf.profiling import (
            profiling_enabled,
            record_sdf_polygon,
        )

        if profiling_enabled():
            record_sdf_polygon(m, elapsed_s=elapsed)
    except Exception:
        pass
    return out
