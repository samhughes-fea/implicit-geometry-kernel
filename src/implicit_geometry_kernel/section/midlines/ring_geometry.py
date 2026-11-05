"""Closed-ring and open-polyline geometry helpers for shell strip topology."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

_MIN_SEGMENT_LENGTH_M: float = 1e-12


def closed_ring_edge_data(
    verts: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """Per-edge length, arc offset at each vertex start, and perimeter ``L``."""
    v = np.asarray(verts, dtype=float)
    nxt = np.roll(v, -1, axis=0)
    edge_len = np.linalg.norm(nxt - v, axis=1)
    L = float(np.sum(edge_len))
    edge_start = np.zeros(int(v.shape[0]), dtype=float)
    if v.shape[0] > 1:
        edge_start[1:] = np.cumsum(edge_len[:-1])
    return edge_len, edge_start, L


def merge_sorted_s_values(values: list[float], tol_m: float) -> list[float]:
    if not values:
        return []
    vals = sorted(float(v) for v in values)
    out = [vals[0]]
    for v in vals[1:]:
        if abs(v - out[-1]) > tol_m:
            out.append(v)
    return out


def nearest_on_closed_polyline(
    pt: NDArray[np.float64],
    verts: NDArray[np.float64],
    edge_len: NDArray[np.float64],
    edge_start: NDArray[np.float64],
    L: float,
) -> tuple[float, NDArray[np.float64]]:
    """Closest point on a closed polygonal ring; returns ``(s, yz)`` with ``s`` in ``[0, L)``."""
    s_arr, yz_arr = nearest_on_closed_polyline_many(
        np.asarray(pt, dtype=float).reshape(1, 2),
        verts,
        edge_len,
        edge_start,
        L,
    )
    return float(s_arr[0]), yz_arr[0].astype(float).copy()


def nearest_on_closed_polyline_many(
    pts: NDArray[np.float64],
    verts: NDArray[np.float64],
    edge_len: NDArray[np.float64],
    edge_start: NDArray[np.float64],
    L: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Batch closest points on a closed ring; returns ``(s[M], yz[M,2])``."""
    p = np.asarray(pts, dtype=float).reshape(-1, 2)
    v = np.asarray(verts, dtype=float)
    if p.shape[0] == 0:
        return np.zeros(0, dtype=float), np.zeros((0, 2), dtype=float)
    if v.shape[0] < 1:
        raise ValueError("closed polyline requires at least one vertex")
    a = v[None, :, :]
    b = np.roll(v, -1, axis=0)[None, :, :]
    e = b - a
    elen2 = np.sum(e * e, axis=2) + 1e-30
    t = np.clip(np.sum((p[:, None, :] - a) * e, axis=2) / elen2, 0.0, 1.0)
    q = a + t[:, :, None] * e
    d2 = np.sum((p[:, None, :] - q) ** 2, axis=2)
    i_best = np.argmin(d2, axis=1)
    m = np.arange(p.shape[0])
    yz = q[m, i_best]
    s = np.asarray(edge_start, dtype=float)[i_best] + t[m, i_best] * np.asarray(
        edge_len, dtype=float
    )[i_best]
    if L > 1e-15:
        s = np.mod(s, L)
    return s.astype(float), yz.astype(float)


def yz_at_s_on_closed_ring(
    s: float,
    verts: NDArray[np.float64],
    edge_len: NDArray[np.float64],
    edge_start: NDArray[np.float64],
    L: float,
) -> NDArray[np.float64]:
    if L <= 1e-15:
        return verts[0].astype(float).copy()
    s_mod = float(s % L)
    n = int(verts.shape[0])
    for i in range(n):
        lo = float(edge_start[i])
        hi = lo + float(edge_len[i])
        if lo <= s_mod <= hi + 1e-12 or (i == n - 1 and s_mod >= lo - 1e-12):
            el = float(edge_len[i])
            t = 0.0 if el < 1e-30 else float(np.clip((s_mod - lo) / el, 0.0, 1.0))
            j = (i + 1) % n
            return verts[i] + t * (verts[j] - verts[i])
    return verts[0].astype(float).copy()


def arc_polyline_no_wrap(
    s_lo: float,
    s_hi: float,
    verts: NDArray[np.float64],
    edge_len: NDArray[np.float64],
    edge_start: NDArray[np.float64],
    L: float,
) -> NDArray[np.float64]:
    """Polyline samples from ``s_lo`` to ``s_hi`` along the ring (``s_lo < s_hi``)."""
    n = int(verts.shape[0])
    out_pts: list[NDArray[np.float64]] = []
    out_pts.append(yz_at_s_on_closed_ring(s_lo, verts, edge_len, edge_start, L))
    for j in range(1, n):
        sj = float(edge_start[j])
        if s_lo + 1e-9 < sj < s_hi - 1e-9:
            out_pts.append(np.asarray(verts[j], dtype=float).copy())
    out_pts.append(yz_at_s_on_closed_ring(s_hi, verts, edge_len, edge_start, L))
    return np.stack(out_pts, axis=0)


def arc_polyline_wrap(
    s_lo: float,
    s_hi: float,
    verts: NDArray[np.float64],
    edge_len: NDArray[np.float64],
    edge_start: NDArray[np.float64],
    L: float,
) -> NDArray[np.float64]:
    """Arc from ``s_lo`` forward through ``L`` and ``0`` to ``s_hi`` (``s_lo > s_hi`` on the circle)."""
    a = arc_polyline_no_wrap(s_lo, L, verts, edge_len, edge_start, L)
    b = arc_polyline_no_wrap(0.0, s_hi, verts, edge_len, edge_start, L)
    if len(a) and len(b) and np.allclose(a[-1], b[0], atol=1e-9, rtol=0.0):
        return np.vstack([a[:-1], b])
    return np.vstack([a, b])


def open_polyline_prefix_lengths(pts: NDArray[np.float64]) -> NDArray[np.float64]:
    """Cumulative distance along open polyline ``pts``; ``prefix[i]`` = arc to vertex ``i``."""
    n = int(pts.shape[0])
    if n < 2:
        return np.zeros(max(n, 0), dtype=float)
    ds = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    prefix = np.zeros(n, dtype=float)
    prefix[1:] = np.cumsum(ds)
    return prefix


def intersect_open_segments_2d(
    a0: NDArray[np.float64],
    a1: NDArray[np.float64],
    b0: NDArray[np.float64],
    b1: NDArray[np.float64],
    *,
    eps: float = 1e-9,
) -> tuple[float, float, NDArray[np.float64]] | None:
    """Strict-interior intersection of segment ``a0a1`` with ``b0b1``; ``ta,tb`` in ``(0,1)``."""
    r = (a1.astype(float) - a0).reshape(2)
    svec = (b1.astype(float) - b0).reshape(2)
    rxs = float(r[0] * svec[1] - r[1] * svec[0])
    denom_tol = 1e-15 * max(
        1.0,
        float(np.linalg.norm(r)) * float(np.linalg.norm(svec)),
    )
    if abs(rxs) < denom_tol:
        return None
    qp = (b0.astype(float) - a0.astype(float)).reshape(2)
    ta = float((qp[0] * svec[1] - qp[1] * svec[0]) / rxs)
    tb = float((qp[0] * r[1] - qp[1] * r[0]) / rxs)
    if ta <= eps or ta >= 1.0 - eps or tb <= eps or tb >= 1.0 - eps:
        return None
    pt = (a0.astype(float) + ta * r).astype(float)
    return ta, tb, pt


def intersect_open_polylines_2d(
    a_pts: NDArray[np.float64],
    b_pts: NDArray[np.float64],
    *,
    eps: float = 1e-9,
) -> list[tuple[int, int, float, float, NDArray[np.float64]]]:
    """All strict-interior edge×edge intersections of two open polylines.

    Returns ``(i_edge_a, j_edge_b, ta, tb, point)`` with ``ta,tb`` in
    ``(eps, 1-eps)``, matching :func:`intersect_open_segments_2d` semantics.

    Both edge indices are reported because a junction has two participants: a
    caller that has to place a shared node on *both* polylines needs to know
    which edge of each one the single solved crossing lands on. Recovering the
    second index afterwards would mean searching for the point with a tolerance,
    which is the failure mode this package keeps designing out.
    """
    a = np.asarray(a_pts, dtype=float)
    b = np.asarray(b_pts, dtype=float)
    if a.ndim != 2 or a.shape[0] < 2 or b.ndim != 2 or b.shape[0] < 2:
        return []
    a0 = a[:-1]
    a1 = a[1:]
    b0 = b[:-1]
    b1 = b[1:]
    r = a1 - a0  # (Na, 2)
    svec = b1 - b0  # (Nb, 2)
    r_norm = np.linalg.norm(r, axis=1)
    s_norm = np.linalg.norm(svec, axis=1)
    rxs = r[:, None, 0] * svec[None, :, 1] - r[:, None, 1] * svec[None, :, 0]
    denom_tol = 1e-15 * np.maximum(1.0, r_norm[:, None] * s_norm[None, :])
    parallel = np.abs(rxs) < denom_tol
    qp = b0[None, :, :] - a0[:, None, :]
    # Avoid divide-by-zero; parallel rows are discarded below.
    rxs_safe = np.where(parallel, 1.0, rxs)
    ta = (qp[..., 0] * svec[None, :, 1] - qp[..., 1] * svec[None, :, 0]) / rxs_safe
    tb = (qp[..., 0] * r[:, None, 1] - qp[..., 1] * r[:, None, 0]) / rxs_safe
    ok = (
        (~parallel)
        & (ta > eps)
        & (ta < 1.0 - eps)
        & (tb > eps)
        & (tb < 1.0 - eps)
    )
    ii, jj = np.nonzero(ok)
    out: list[tuple[int, int, float, float, NDArray[np.float64]]] = []
    for i, j in zip(ii.tolist(), jj.tolist()):
        ta_i = float(ta[i, j])
        tb_j = float(tb[i, j])
        pt = (a0[i] + ta_i * r[i]).astype(float)
        out.append((int(i), int(j), ta_i, tb_j, pt.reshape(2)))
    return out


def intersect_web_axis_open_polyline_2d(
    web_pts: NDArray[np.float64],
    cap_pts: NDArray[np.float64],
    *,
    endpoint_tol_m: float = 1e-9,
) -> list[tuple[int, int, float, float, NDArray[np.float64]]]:
    """Interior hits of the web axis (``web[0]→web[-1]``) against an open polyline.

    Returns ``(i_cap_edge, j_web_edge, ta, tb, point)`` in the same shape as
    :func:`intersect_open_polylines_2d`. Cap **endpoints** are excluded (those are
    cap-end junctions closed at emit). ``ta`` is the parameter on the cap edge;
    ``tb`` is the parameter on the discrete web edge that carries the same point
    along the axis (for insertion into the sampled web polyline).

    This is the production interior-junction solver (SG-009). Prefer it over
    :func:`intersect_open_polylines_2d` for topology truth: polyline×polyline
    detection cannot verify its own work once a shared node has been inserted
    (tangency trap, SG-005).
    """
    web = np.asarray(web_pts, dtype=float)
    cap = np.asarray(cap_pts, dtype=float)
    if web.ndim != 2 or web.shape[0] < 2 or cap.ndim != 2 or cap.shape[0] < 2:
        return []

    origin = web[0]
    span = web[-1] - web[0]
    length = float(np.linalg.norm(span))
    if length < 1e-12:
        return []
    direction = span / length

    seg_start = cap[:-1]
    seg = cap[1:] - seg_start
    denominator = seg[:, 0] * direction[1] - seg[:, 1] * direction[0]
    usable = np.abs(denominator) > 1e-14
    safe = np.where(usable, denominator, 1.0)
    rel = origin - seg_start
    u = np.where(usable, (rel[:, 0] * direction[1] - rel[:, 1] * direction[0]) / safe, -9.0)
    t_axis = np.where(
        usable, (rel[:, 0] * seg[:, 1] - rel[:, 1] * seg[:, 0]) / safe, -9.0
    )
    on_both = (
        usable
        & (u >= -1e-9)
        & (u <= 1.0 + 1e-9)
        & (t_axis >= -1e-9)
        & (t_axis <= length + 1e-9)
    )

    # Axis parameter of each web vertex, for mapping hits onto discrete edges.
    t_vert = np.sum((web - origin) * direction, axis=1)

    out: list[tuple[int, int, float, float, NDArray[np.float64]]] = []
    end_tol = float(endpoint_tol_m)
    for i in np.flatnonzero(on_both):
        ta = float(np.clip(u[i], 0.0, 1.0))
        t_hit = float(np.clip(t_axis[i], 0.0, length))
        point = (seg_start[i] + ta * seg[i]).astype(float).reshape(2)
        at_cap_end = min(
            float(np.linalg.norm(point - cap[0])),
            float(np.linalg.norm(point - cap[-1])),
        )
        if at_cap_end < end_tol:
            continue

        j_web = -1
        tb = 0.0
        for j in range(int(web.shape[0]) - 1):
            t0 = float(t_vert[j])
            t1 = float(t_vert[j + 1])
            span_e = t1 - t0
            if abs(span_e) < 1e-15:
                if abs(t_hit - t0) <= 1e-9:
                    j_web = j
                    tb = 0.0
                    break
                continue
            # Inclusive so a hit at an existing web vertex still maps.
            lo, hi = (t0, t1) if span_e > 0.0 else (t1, t0)
            if t_hit < lo - 1e-9 or t_hit > hi + 1e-9:
                continue
            tb = (t_hit - t0) / span_e
            tb = float(np.clip(tb, 0.0, 1.0))
            j_web = j
            break
        if j_web < 0:
            continue
        out.append((int(i), int(j_web), ta, tb, point))
    return out


def open_polyline_min_segment_length(pts: NDArray[np.float64]) -> float:
    if pts.shape[0] < 2:
        return float("inf")
    return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).min())


def point_to_open_polyline_nearest(
    p: NDArray[np.float64], w: NDArray[np.float64]
) -> tuple[float, NDArray[np.float64]]:
    """Nearest distance and point from 2-D point ``p`` to open polyline ``w``."""
    q = np.asarray(p, dtype=float).reshape(2)
    pts = np.asarray(w, dtype=float)
    if pts.shape[0] < 2:
        return float("inf"), q.copy()
    a = pts[:-1]
    b = pts[1:]
    ab = b - a
    l2 = np.sum(ab * ab, axis=1) + 1e-30
    t = np.clip(np.sum((q - a) * ab, axis=1) / l2, 0.0, 1.0)
    c = a + t[:, None] * ab
    d = np.linalg.norm(q - c, axis=1)
    j = int(np.argmin(d))
    return float(d[j]), np.asarray(c[j], dtype=float).reshape(2)


def normal_2d_ccw(t: NDArray[np.float64]) -> NDArray[np.float64]:
    """Unit normal: 90° CCW rotation of a 2-D tangent (same convention as MITC4 assembly)."""
    tx, ty = float(t[0]), float(t[1])
    n = np.array([-ty, tx], dtype=float)
    ln = float(np.linalg.norm(n))
    if ln < 1e-30:
        return np.array([0.0, 1.0], dtype=float)
    return n / ln


__all__ = [
    "_MIN_SEGMENT_LENGTH_M",
    "closed_ring_edge_data",
    "merge_sorted_s_values",
    "nearest_on_closed_polyline",
    "nearest_on_closed_polyline_many",
    "yz_at_s_on_closed_ring",
    "arc_polyline_no_wrap",
    "arc_polyline_wrap",
    "open_polyline_prefix_lengths",
    "intersect_open_segments_2d",
    "intersect_open_polylines_2d",
    "intersect_web_axis_open_polyline_2d",
    "open_polyline_min_segment_length",
    "point_to_open_polyline_nearest",
    "normal_2d_ccw",
]
