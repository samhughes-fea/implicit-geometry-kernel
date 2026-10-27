"""Resample cap interiors via inward rays from the skin master ring."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from numpy.typing import NDArray

from blade_geometry_model.section_contract.shell_midline_strip import ShellMidlineStrip

from .ring_geometry import (
    _MIN_SEGMENT_LENGTH_M,
    normal_2d_ccw,
    open_polyline_prefix_lengths,
)

_SKIN_RAY_T_MIN_M: float = 1e-9


def _surface_side(value: object) -> str | None:
    """Normalise a cap ``surface`` or skin ``skin_side`` to suction/pressure."""
    text = str(getattr(value, "value", value) or "").lower()
    if "suction" in text or text in ("ss", "upper"):
        return "suction_surface"
    if "pressure" in text or text in ("ps", "lower"):
        return "pressure_surface"
    return None


def concat_skin_ring_vertices(
    strips: tuple[ShellMidlineStrip, ...],
    *,
    with_sides: bool = False,
):
    """Concatenate all skin ``midline_twisted_pitch_axis_blade_frame`` strips
    in list order; dedupe shared joints.

    With ``with_sides`` the return is ``(ring, sides)`` where ``sides`` labels
    each ring vertex suction/pressure, so a cap can be resampled from the skin
    it actually lies under rather than from whichever vertices happen to fall
    inside a bounding box.
    """
    chunks: list[NDArray[np.float64]] = []
    side_chunks: list[NDArray[np.object_]] = []
    for st in strips:
        if st.kind != "skin":
            continue
        arr = np.asarray(st.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if arr.shape[0] < 2:
            continue
        if chunks and np.allclose(chunks[-1][-1], arr[0], atol=1e-9, rtol=0.0):
            arr = arr[1:]
        if arr.shape[0] == 0:
            continue
        chunks.append(arr)
        side_chunks.append(
            np.full(arr.shape[0], _surface_side(st.skin_side), dtype=object)
        )
    if not chunks:
        empty = np.empty((0, 2), dtype=float)
        return (empty, np.empty(0, dtype=object)) if with_sides else empty
    ring = np.vstack(chunks)
    if not with_sides:
        return ring
    return ring, np.concatenate(side_chunks)


def skin_ring_vertex_tangents(ring: NDArray[np.float64]) -> NDArray[np.float64]:
    """Unit tangent at each ring vertex; ring treated closed (last → first edge)."""
    r = np.asarray(ring, dtype=float)
    m = int(r.shape[0])
    if m < 2:
        return np.zeros((0, 2), dtype=float)
    prev_v = np.roll(r, 1, axis=0)
    next_v = np.roll(r, -1, axis=0)
    d = next_v - prev_v
    ln = np.linalg.norm(d, axis=1)
    out = np.zeros((m, 2), dtype=float)
    ok = ln >= 1e-30
    out[ok] = d[ok] / ln[ok, None]
    out[~ok] = np.array([1.0, 0.0], dtype=float)
    return out


def _ray_segment_hit_first(
    p: NDArray[np.float64],
    v: NDArray[np.float64],
    a: NDArray[np.float64],
    b: NDArray[np.float64],
    *,
    t_min: float = _SKIN_RAY_T_MIN_M,
) -> tuple[float, NDArray[np.float64]] | None:
    """First hit of ray ``p + t v`` (``t >= t_min``) on segment ``ab``; ``v`` unit."""
    w = b - a
    mat = np.column_stack([v, -w])
    det = float(np.linalg.det(mat))
    if abs(det) < 1e-18:
        return None
    rhs = a - p
    tu = np.linalg.solve(mat, rhs)
    t = float(tu[0])
    u = float(tu[1])
    if t < t_min or u < -1e-10 or u > 1.0 + 1e-10:
        return None
    uu = float(np.clip(u, 0.0, 1.0))
    hit = a + uu * w
    return t, hit


def _ray_hits_open_polyline(
    p: NDArray[np.float64],
    v: NDArray[np.float64],
    cap_pts: NDArray[np.float64],
) -> list[tuple[float, float, NDArray[np.float64]]]:
    """All ray hits on open ``cap_pts``; each entry ``(t_ray, s_cap, hit)`` with smallest ``t`` first."""
    pts = np.asarray(cap_pts, dtype=float)
    n = int(pts.shape[0])
    if n < 2:
        return []
    p2 = np.asarray(p, dtype=float).reshape(2)
    v2 = np.asarray(v, dtype=float).reshape(2)
    a = pts[:-1]
    b = pts[1:]
    w = b - a
    vx, vy = float(v2[0]), float(v2[1])
    # det([[vx,-wx],[vy,-wy]]) = -vx*wy + wx*vy
    det = -vx * w[:, 1] + w[:, 0] * vy
    rhs = a - p2
    # Avoid divide-by-zero; mark parallel edges invalid.
    ok_det = np.abs(det) >= 1e-18
    det_safe = np.where(ok_det, det, 1.0)
    t = (-rhs[:, 0] * w[:, 1] + w[:, 0] * rhs[:, 1]) / det_safe
    u = (vx * rhs[:, 1] - vy * rhs[:, 0]) / det_safe
    t_min = _SKIN_RAY_T_MIN_M
    ok = ok_det & (t >= t_min) & (u >= -1e-10) & (u <= 1.0 + 1e-10)
    if not np.any(ok):
        return []
    uu = np.clip(u[ok], 0.0, 1.0)
    hit = a[ok] + uu[:, None] * w[ok]
    prefix = open_polyline_prefix_lengths(pts)
    u_len = np.linalg.norm(w[ok], axis=1)
    s_cap = prefix[:-1][ok] + np.where(
        u_len < 1e-30, 0.0, uu * u_len
    )
    order = np.argsort(t[ok])
    hits: list[tuple[float, float, NDArray[np.float64]]] = []
    t_ok = t[ok]
    for idx in order:
        hits.append(
            (float(t_ok[idx]), float(s_cap[idx]), hit[idx].astype(float).copy())
        )
    return hits


def pick_inward_ray_direction(
    p_skin: NDArray[np.float64],
    n_ccw: NDArray[np.float64],
    cap_centroid: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Choose ``+/- n_ccw`` so the ray points generally toward the cap interior."""
    n1 = np.asarray(n_ccw, dtype=float)
    n2 = -n1
    if float(np.dot(n1, cap_centroid - p_skin)) >= float(np.dot(n2, cap_centroid - p_skin)):
        return n1
    return n2


def _resample_one_cap_interior_from_skin_rays(
    cap_strip: ShellMidlineStrip,
    skin_ring: NDArray[np.float64],
    skin_tangents: NDArray[np.float64],
    *,
    skin_sides: NDArray[np.object_] | None = None,
    cap_side: str | None = None,
    junction_seeds: list[NDArray[np.float64]] | None = None,
) -> ShellMidlineStrip | None:
    """Return replaced cap strip, or ``None`` to keep original.

    ``junction_seeds`` are points that must survive the rebuild (typically
    pre-resample web-axis × cap crossings). They are re-inserted onto the
    resampled polyline by nearest-edge projection so later topology can
    find conformal T-junctions without rediscovering tangencies.
    """
    pts = np.asarray(cap_strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
    if pts.shape[0] < 3:
        return None

    p0 = pts[0].copy()
    p1 = pts[-1].copy()
    c_mid = np.asarray(np.mean(pts, axis=0), dtype=float).reshape(2)

    prefix_full = open_polyline_prefix_lengths(pts)
    arc_total = float(prefix_full[-1])
    if arc_total < 1e-9:
        return None
    s_eps = max(1e-6 * arc_total, 1e-7)
    # Target interior spacing follows the incoming cap polyline resolution so
    # skin-ray hits do not over-cluster and create pathological short segments.
    target_step = arc_total / max(float(pts.shape[0] - 1), 1.0)

    collected: list[tuple[float, NDArray[np.float64]]] = []

    m = int(skin_ring.shape[0])
    for i in range(m):
        # Which skin a cap belongs under is declared on the strip, so there is
        # nothing to infer geometrically. Vertices from the other surface are
        # excluded because a ray aimed at the cap centroid from the far side can
        # still strike it; every remaining hit is filtered by arc position below,
        # which is what actually bounds the resample.
        if skin_sides is not None and skin_sides[i] != cap_side:
            continue
        q = skin_ring[i]
        ti = skin_tangents[i]
        n_ccw = normal_2d_ccw(ti)
        v = pick_inward_ray_direction(q, n_ccw, c_mid)
        hits = _ray_hits_open_polyline(q, v, pts)
        for t_ray, s_cap, hit in hits:
            if s_cap <= s_eps or s_cap >= arc_total - s_eps:
                continue
            collected.append((s_cap, hit))
            break

    if len(collected) < 1 and not junction_seeds:
        return None

    collected.sort(key=lambda x: x[0])
    merged: list[tuple[float, NDArray[np.float64]]] = []
    merge_tol = max(0.35 * target_step, 1e-7)
    for s_cap, hit in collected:
        if merged and abs(s_cap - merged[-1][0]) <= merge_tol:
            continue
        merged.append((s_cap, hit))

    seq: list[NDArray[np.float64]] = [p0.copy()]
    min_keep = max(0.5 * target_step, _MIN_SEGMENT_LENGTH_M)
    for _s_cap, hit in merged:
        if float(np.linalg.norm(hit - seq[-1])) >= min_keep:
            seq.append(hit.copy())
    if float(np.linalg.norm(p1 - seq[-1])) >= min_keep:
        seq.append(p1.copy())
    else:
        seq[-1] = p1.copy()

    new_pts = np.stack(seq, axis=0)
    if new_pts.shape[0] < 2:
        return None
    if not np.allclose(new_pts[0], p0, atol=1e-12, rtol=0.0) or not np.allclose(
        new_pts[-1], p1, atol=1e-12, rtol=0.0
    ):
        return None

    if junction_seeds:
        new_pts = _reinsert_junction_seeds(new_pts, junction_seeds, min_seg=min_keep)

    return replace(cap_strip, midline_twisted_pitch_axis_blade_frame=new_pts)


def _reinsert_junction_seeds(
    pts: NDArray[np.float64],
    seeds: list[NDArray[np.float64]],
    *,
    min_seg: float,
) -> NDArray[np.float64]:
    """Insert seed points that already lie near ``pts`` (no relocation)."""
    out = np.asarray(pts, dtype=float)
    for seed in seeds:
        p = np.asarray(seed, dtype=float).reshape(2)
        if float(np.linalg.norm(out - p, axis=1).min()) <= 1e-9:
            continue
        # Project seed onto nearest edge; if far, skip (seed was for old shape).
        a = out[:-1]
        b = out[1:]
        ab = b - a
        l2 = np.sum(ab * ab, axis=1) + 1e-30
        t = np.sum((p[None, :] - a) * ab, axis=1) / l2
        best_i = -1
        best_d = float("inf")
        best_t = 0.0
        for i in range(len(a)):
            ti = float(np.clip(t[i], 0.0, 1.0))
            c = a[i] + ab[i] * ti
            d = float(np.linalg.norm(p - c))
            if d < best_d:
                best_d = d
                best_i = i
                best_t = ti
        if best_i < 0 or best_d > max(1e-3, 5.0 * min_seg):
            continue
        # Prefer the exact seed when already on the polyline within a hair.
        insert_p = p if best_d <= 1e-9 else (a[best_i] + ab[best_i] * best_t)
        rows: list[NDArray[np.float64]] = []
        for i in range(out.shape[0] - 1):
            rows.append(out[i].copy())
            if i == best_i:
                if float(np.linalg.norm(insert_p - rows[-1])) >= min_seg and float(
                    np.linalg.norm(out[i + 1] - insert_p)
                ) >= min_seg:
                    rows.append(np.asarray(insert_p, dtype=float).reshape(2))
        rows.append(out[-1].copy())
        out = np.stack(rows, axis=0)
    return out


def _collect_cap_web_junction_seeds(
    strips: tuple[ShellMidlineStrip, ...],
) -> dict[int, list[NDArray[np.float64]]]:
    """Pre-resample interior crossings via junction_primitives (axis × cap)."""
    from blade_geometry_model.section_geometry.midlines.junction_primitives import (
        intersect_axis_cap_polyline,
    )
    from blade_geometry_model.section_geometry.midlines.ring_geometry import (
        open_polyline_prefix_lengths,
    )

    webs = [
        np.asarray(st.midline_twisted_pitch_axis_blade_frame, dtype=float)
        for st in strips
        if st.kind == "web"
        and np.asarray(st.midline_twisted_pitch_axis_blade_frame).shape[0] >= 2
    ]
    seeds: dict[int, list[NDArray[np.float64]]] = {}
    for idx, st in enumerate(strips):
        if st.kind != "cap":
            continue
        pts = np.asarray(st.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if pts.shape[0] < 2:
            continue
        prefix = open_polyline_prefix_lengths(pts)
        arc = float(prefix[-1])
        merge = 1e-7 * max(arc, 1e-15)
        collected: list[NDArray[np.float64]] = []
        for w in webs:
            for i, _j, ta, _tb, p_hit in intersect_axis_cap_polyline(w, pts):
                seg = float(np.linalg.norm(pts[i + 1] - pts[i]))
                s_hit = float(prefix[i] + ta * seg)
                if s_hit <= merge or s_hit >= arc - merge:
                    continue
                collected.append(np.asarray(p_hit, dtype=float).reshape(2))
        if collected:
            seeds[idx] = collected
    return seeds


def resample_cap_interiors_from_skin_rays(
    strips: tuple[ShellMidlineStrip, ...],
) -> tuple[ShellMidlineStrip, ...]:
    """Replace interior vertices of each cap using inward rays from skin master ring.

    Cap **endpoints** are left unchanged — they are already at their true extent
    on emit, and nothing downstream relocates them. Pre-resample web×cap
    **junction seeds** are collected and re-inserted after the rebuild so later
    topology does not have to rediscover crossings that resampling would erase.

    If resampling yields too few hits for a cap, that strip is left unchanged.

    Which skin a cap sits under is taken from the strips' declared
    ``surface`` / ``skin_side``, not from a bounding box. The bounding-box
    pre-filter this replaced carried a fixed margin that could be narrower than
    the skin-to-cap separation ``(skin_t + cap_t)/2``, in which case it discarded
    nearly every skin vertex above the cap and collapsed its interior to a chord.
    """
    skin_ring, skin_sides = concat_skin_ring_vertices(strips, with_sides=True)
    if skin_ring.shape[0] < 3:
        return strips

    skin_tangents = skin_ring_vertex_tangents(skin_ring)
    if skin_tangents.shape[0] != skin_ring.shape[0]:
        return strips

    seeds_by_cap = _collect_cap_web_junction_seeds(strips)

    out: list[ShellMidlineStrip] = []
    for idx, st in enumerate(strips):
        if st.kind != "cap":
            out.append(st)
            continue
        new_st = _resample_one_cap_interior_from_skin_rays(
            st,
            skin_ring,
            skin_tangents,
            skin_sides=skin_sides,
            cap_side=_surface_side(st.surface),
            junction_seeds=seeds_by_cap.get(idx),
        )
        out.append(new_st if new_st is not None else st)
    return tuple(out)


__all__ = [
    "concat_skin_ring_vertices",
    "skin_ring_vertex_tangents",
    "pick_inward_ray_direction",
    "resample_cap_interiors_from_skin_rays",
]
