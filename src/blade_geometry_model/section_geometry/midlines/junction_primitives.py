"""Analytic junction primitives (axis × offset locus / final cap polyline).

Production web feet use :func:`solve_web_foot_on_skin_mid` (same
``offset_locus`` join contract as constructive polygons). Cap interiors close
post-resample via :func:`intersect_axis_cap_polyline` (web axis × final cap —
the discrete offset locus under SG-008). When an untwisted airfoil and cap mid
depth are supplied, :func:`solve_interior_cap_web_hits` prefers
axis × ``locus_at_depth`` hits that land on the final cap polyline.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from blade_geometry_model.section_geometry.midlines.offset_ring import (
    locus_at_depth,
)
from blade_geometry_model.section_geometry.midlines.ring_geometry import (
    intersect_web_axis_open_polyline_2d,
)


def intersect_axis_offset_ring(
    axis_xy: NDArray[np.float64] | list,
    airfoil,
    depth: float,
) -> list[NDArray[np.float64]]:
    """Intersect a 2-point axis with the offset locus at ``depth``.

    Returns hit points (may be empty). The ring is closed for edge coverage.
    """
    axis = np.asarray(axis_xy, dtype=float)
    if axis.ndim != 2 or axis.shape[0] < 2 or axis.shape[1] != 2:
        raise ValueError(f"axis_xy must be (N,2) with N>=2; got {axis.shape}")
    ring = np.asarray(locus_at_depth(airfoil, float(depth)), dtype=float)
    if ring.shape[0] < 2:
        return []
    if float(np.hypot(*(ring[0] - ring[-1]))) > 1e-15:
        ring = np.vstack([ring, ring[:1]])
    hits = intersect_web_axis_open_polyline_2d(axis, ring, endpoint_tol_m=0.0)
    return [np.asarray(h[-1], dtype=float) for h in hits]


def solve_web_foot_on_skin_mid(
    axis_p0,
    axis_p1,
    airfoil,
    skin_thickness: float,
) -> list[NDArray[np.float64]]:
    """Web-foot candidates: axis × skin mid-surface (offset depth = t/2)."""
    axis = np.asarray([axis_p0, axis_p1], dtype=float)
    return intersect_axis_offset_ring(axis, airfoil, 0.5 * float(skin_thickness))


def intersect_axis_cap_polyline(
    web_pts: NDArray[np.float64],
    cap_pts: NDArray[np.float64],
    *,
    endpoint_tol_m: float = 1e-9,
) -> list[tuple[int, int, float, float, NDArray[np.float64]]]:
    """Production interior detector: web axis × final cap polyline (SG-009).

    The final cap polyline is the discrete mid-surface locus (SG-008); this is
    the same geometric contract as :func:`intersect_axis_offset_ring`, evaluated
    on the band-clipped sample rather than the full airfoil ring.
    """
    return intersect_web_axis_open_polyline_2d(
        web_pts, cap_pts, endpoint_tol_m=endpoint_tol_m
    )


def solve_interior_cap_web_hits(
    web_pts: NDArray[np.float64],
    cap_pts: NDArray[np.float64],
    *,
    airfoil=None,
    depth: float | None = None,
    on_cap_tol_m: float = 1e-6,
    endpoint_tol_m: float = 1e-9,
) -> list[tuple[int, int, float, float, NDArray[np.float64]]]:
    """Interior crossings for topology.

    Prefer axis × offset locus at ``depth`` when ``airfoil`` is given (same-frame
    as the strips — typically twist ≈ 0). Hits must land on the final cap
    polyline within ``on_cap_tol_m``. Otherwise fall back to
    :func:`intersect_axis_cap_polyline`.
    """
    poly_hits = intersect_axis_cap_polyline(
        web_pts, cap_pts, endpoint_tol_m=endpoint_tol_m
    )
    if airfoil is None or depth is None:
        return poly_hits

    axis = np.asarray(web_pts, dtype=float)
    ring_hits = intersect_axis_offset_ring(axis, airfoil, float(depth))
    if not ring_hits:
        return poly_hits

    # Map each ring hit onto the final cap polyline; keep only on-cap landings.
    pts = np.asarray(cap_pts, dtype=float)
    if pts.shape[0] < 2:
        return poly_hits
    a0 = pts[:-1]
    a1 = pts[1:]
    ab = a1 - a0
    l2 = np.sum(ab * ab, axis=1) + 1e-30
    out: list[tuple[int, int, float, float, NDArray[np.float64]]] = []
    for p_ring in ring_hits:
        q = np.asarray(p_ring, dtype=float).reshape(2)
        t_raw = np.sum((q[None, :] - a0) * ab, axis=1) / l2
        best = None
        best_gap = float("inf")
        for i in range(int(pts.shape[0]) - 1):
            ti = float(np.clip(t_raw[i], 0.0, 1.0))
            c = a0[i] + ab[i] * ti
            gap = float(np.linalg.norm(q - c))
            if gap < best_gap:
                best_gap = gap
                best = (i, ti, c)
        if best is None or best_gap > float(on_cap_tol_m):
            continue
        i_cap, t_cap, p_on = best
        # Match return order of intersect_web_axis_open_polyline_2d:
        # (cap_edge, web_edge, t_cap, t_web, point).
        w_hit = None
        for i_c, j_w, ta, tb, p_h in poly_hits:
            if i_c == i_cap and abs(float(ta) - t_cap) < 1e-6:
                w_hit = (
                    i_c,
                    j_w,
                    float(ta),
                    float(tb),
                    np.asarray(p_on, dtype=float),
                )
                break
        if w_hit is None:
            w0, w1 = axis[0], axis[-1]
            w_ab = w1 - w0
            w_l2 = float(np.dot(w_ab, w_ab)) + 1e-30
            t_w = float(np.clip(np.dot(p_on - w0, w_ab) / w_l2, 0.0, 1.0))
            w_hit = (i_cap, 0, float(t_cap), t_w, np.asarray(p_on, dtype=float))
        out.append(w_hit)

    return out if out else poly_hits


__all__ = [
    "intersect_axis_cap_polyline",
    "intersect_axis_offset_ring",
    "solve_interior_cap_web_hits",
    "solve_web_foot_on_skin_mid",
]
