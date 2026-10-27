"""Insert cap midline vertices at web T-junctions."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from numpy.typing import NDArray

from blade_geometry_model.section_contract.shell_midline_strip import ShellMidlineStrip

from ..midlines.junction_primitives import solve_interior_cap_web_hits
from ..midlines.ring_geometry import (
    _MIN_SEGMENT_LENGTH_M,
    open_polyline_min_segment_length,
    open_polyline_prefix_lengths,
)

_CAP_WEB_JUNCTION_MERGE_ARC_FRAC: float = 1e-7


def _insert_points_on_segments(
    pts: NDArray[np.float64],
    insertions: list[tuple[int, float, NDArray[np.float64]]],
    *,
    vertex_coincide_tol_m: float,
    min_seg: float,
) -> NDArray[np.float64] | None:
    """Add vertices that already lie on ``pts``, or ``None`` if nothing changes.

    Each insertion is ``(edge_index, t, point)`` with ``point`` on edge
    ``edge_index``. Nothing is projected, snapped or moved: the polyline's shape
    is identical afterwards and only its vertex set grows. That exactness is the
    whole point — the same solved crossing goes into both participants of a
    junction, so neither has to find the other by tolerance.
    """
    by_edge: dict[int, list[tuple[float, NDArray[np.float64]]]] = {}
    for edge, t, p in insertions:
        by_edge.setdefault(int(edge), []).append(
            (float(t), np.asarray(p, dtype=float).reshape(2))
        )
    if not by_edge:
        return None

    rows: list[NDArray[np.float64]] = []
    changed = False
    n = int(pts.shape[0])
    for i in range(n - 1):
        rows.append(pts[i].copy())
        for _t, p in sorted(by_edge.get(i, []), key=lambda z: z[0]):
            if float(np.linalg.norm(pts - p, axis=1).min()) <= vertex_coincide_tol_m:
                continue
            if float(np.linalg.norm(p - rows[-1])) < min_seg:
                continue
            if float(np.linalg.norm(pts[i + 1] - p)) < min_seg:
                continue
            rows.append(p.copy())
            changed = True
    rows.append(pts[-1].copy())

    if not changed:
        return None
    new_pts = np.stack(rows, axis=0)
    if new_pts.shape[0] < 2 or open_polyline_min_segment_length(new_pts) < _MIN_SEGMENT_LENGTH_M:
        return None
    return new_pts


def insert_cap_web_interior_junctions(
    strips: tuple[ShellMidlineStrip, ...],
    *,
    vertex_coincide_tol_m: float = 1e-9,
    merge_arc_m: float | None = None,
    min_seg_m: float | None = None,
    foot_proximity_tol_m: float = 0.0,
    airfoil=None,
    cap_mid_depth_m: float | None = None,
) -> tuple[ShellMidlineStrip, ...]:
    """Give cap **and** web a shared node at every interior cap-web T-junction.

    (1) Solve each interior crossing via
    :func:`midlines.junction_primitives.solve_interior_cap_web_hits` —
    preferably axis × offset locus at ``cap_mid_depth_m`` when ``airfoil`` is
    supplied (same frame as the strips), otherwise web axis × final cap
    polyline (SG-009 / SG-008 discrete locus). Each hit is solved once and the
    *same* point is inserted into both strips.
    (2) Optional **proximity** to each shear-web **foot** (endpoint on skin): when
    enabled (``foot_proximity_tol_m > 0``) and the foot is within tolerance of the
    cap polyline, insert the foot ``yz`` as a new vertex. Cap-side only — the foot
    is already a web vertex.

    Layout-wise, **Y = C** systems hit every web this way by construction; **Y = B**
    systems often need it only once **X** is large enough that the cap band
    crosses webs between bay ends.

    Intended after :func:`split_skin_at_junctions` and (when enabled) after
    :func:`resample_cap_interiors_from_skin_rays`, so cap endpoints are final.

    **Why both sides here and not at emit.** ``shell_midline_export`` closes the
    cap-*end* junctions on the web at emit, but it cannot close interior ones:
    the cap it would compute against is not the cap that survives resampling.
    Nor can the web side be done first and the cap side detected afterwards —
    placing a web vertex exactly on the cap turns a polyline×polyline crossing
    into a tangency (SG-005 ``4D-F`` swap). Axis×cap detection does not have that
    trap, but insertions are still held and applied afterwards so both sides get
    the same solved point. See ``docs/bug_log.md`` SG-005 / SG-009.
    """
    web_pts: list[NDArray[np.float64]] = []
    web_feet: list[NDArray[np.float64]] = []
    web_strip_index: list[int] = []
    for idx, st in enumerate(strips):
        if st.kind != "web":
            continue
        w = np.asarray(st.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if w.shape[0] >= 2:
            web_pts.append(w)
            web_strip_index.append(idx)
            web_feet.append(np.asarray(w[0], dtype=float).reshape(2))
            web_feet.append(np.asarray(w[-1], dtype=float).reshape(2))

    if not web_pts:
        return strips

    # strip index -> [(web edge index, t along that edge, crossing point)]
    web_insertions: dict[int, list[tuple[int, float, NDArray[np.float64]]]] = {}

    min_seg = float(min_seg_m) if min_seg_m is not None else max(_MIN_SEGMENT_LENGTH_M * 100.0, 1e-9)
    foot_tol = float(foot_proximity_tol_m)
    use_foot_proximity = foot_tol > 0.0

    out: list[ShellMidlineStrip] = []
    for st in strips:
        if st.kind != "cap":
            out.append(st)
            continue
        pts = np.asarray(st.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if pts.shape[0] < 2:
            out.append(st)
            continue
        prefix = open_polyline_prefix_lengths(pts)
        arc_total = float(prefix[-1])
        if arc_total < min_seg * 3.0:
            out.append(st)
            continue
        merge_arc = (
            float(merge_arc_m)
            if merge_arc_m is not None
            else _CAP_WEB_JUNCTION_MERGE_ARC_FRAC * arc_total
        )
        eps_s = max(1e-12 * arc_total, 1e-15)

        hits: list[tuple[float, NDArray[np.float64]]] = []
        # Crossings against this cap, held back until the cap's own polyline is
        # final: a point only earns a web node if the cap really carries it.
        pending_web: list[tuple[int, int, float, NDArray[np.float64]]] = []
        n = int(pts.shape[0])
        seg_lens = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        for w, w_idx in zip(web_pts, web_strip_index, strict=True):
            for i, j, ta, tb, p_hit in solve_interior_cap_web_hits(
                w,
                pts,
                airfoil=airfoil,
                depth=cap_mid_depth_m,
                on_cap_tol_m=max(vertex_coincide_tol_m, 1e-6),
            ):
                seg_len = float(seg_lens[i])
                if seg_len < min_seg:
                    continue
                s_hit = float(prefix[i] + ta * seg_len)
                if s_hit <= merge_arc or s_hit >= arc_total - merge_arc:
                    continue
                p = np.asarray(p_hit, dtype=float).reshape(2)
                hits.append((s_hit, p))
                pending_web.append((w_idx, j, tb, p))

        if use_foot_proximity and web_feet:
            a0 = pts[:-1]
            a1 = pts[1:]
            ab = a1 - a0
            l2 = np.sum(ab * ab, axis=1) + 1e-30
            feet = np.asarray(web_feet, dtype=float).reshape(-1, 2)
            for q in feet:
                t_raw = np.sum((q[None, :] - a0) * ab, axis=1) / l2
                for i in range(n - 1):
                    seg_len = float(seg_lens[i])
                    if seg_len < min_seg:
                        continue
                    eps_line = max(1e-12, min(1e-7, 1e-6 * seg_len / max(seg_len, 1e-15)))
                    ti = float(t_raw[i])
                    if ti <= 0.0:
                        if float(np.linalg.norm(q - a0[i])) > foot_tol:
                            continue
                        t_ins = eps_line
                    elif ti >= 1.0:
                        if float(np.linalg.norm(q - a1[i])) > foot_tol:
                            continue
                        t_ins = 1.0 - eps_line
                    else:
                        c = a0[i] + ab[i] * ti
                        if float(np.linalg.norm(q - c)) > foot_tol:
                            continue
                        t_ins = ti
                    s_hit = float(prefix[i] + t_ins * seg_len)
                    if s_hit <= merge_arc or s_hit >= arc_total - merge_arc:
                        continue
                    if float(np.linalg.norm(pts - q, axis=1).min()) <= vertex_coincide_tol_m:
                        continue
                    hits.append((s_hit, q.copy()))

        # A cap-endpoint-on-web branch used to live here, inserting the nearest
        # web-midline point for any cap end within ``_CAP_ENDPOINT_ON_WEB_TOL_M``
        # of a web. It became unreachable once ``shell_midline_export`` started
        # emitting the web vertex at each cap-end junction: a family-D cap end is
        # then already coincident with a web node, so the ``vertex_coincide_tol_m``
        # guard skipped it, and a family-B/C free end is 38.6-148.6 mm from the
        # nearest web, so the tolerance guard skipped it. Removed after
        # confirming byte-identical output across 108 layout x twist x chord
        # configurations with the branch disabled. See ``docs/bug_log.md`` SG-004.

        def _queue_web_nodes(final_cap_pts: NDArray[np.float64]) -> None:
            """Queue the web half of every crossing the cap actually carries.

            Guards against a half-made junction: if a crossing was dropped from
            the cap (slimmed away as too short, say), putting the node on the
            web alone would move the defect rather than fix it. A crossing that
            was skipped because the cap *already* had a vertex there still
            qualifies — the cap has the node either way.
            """
            for w_idx, j_edge, t_edge, p in pending_web:
                gap = float(np.linalg.norm(final_cap_pts - p, axis=1).min())
                if gap > vertex_coincide_tol_m:
                    continue
                web_insertions.setdefault(w_idx, []).append((j_edge, t_edge, p))

        if not hits:
            _queue_web_nodes(pts)
            out.append(st)
            continue

        hits.sort(key=lambda h: h[0])
        merged: list[tuple[float, NDArray[np.float64]]] = []
        for s_h, p_h in hits:
            if merged and abs(s_h - merged[-1][0]) <= merge_arc:
                continue
            if float(np.linalg.norm(pts - p_h, axis=1).min()) <= vertex_coincide_tol_m:
                continue
            merged.append((s_h, p_h))

        if not merged:
            _queue_web_nodes(pts)
            out.append(st)
            continue

        new_rows: list[NDArray[np.float64]] = []
        mi = 0
        for i in range(n - 1):
            new_rows.append(pts[i].copy())
            s_lo = float(prefix[i])
            s_hi = float(prefix[i + 1])
            while mi < len(merged) and merged[mi][0] < s_hi - eps_s:
                s_h, p_h = merged[mi]
                if s_h > s_lo + eps_s and float(np.linalg.norm(p_h - new_rows[-1])) >= min_seg:
                    new_rows.append(p_h.copy())
                mi += 1
        new_rows.append(pts[-1].copy())

        slim: list[NDArray[np.float64]] = [new_rows[0]]
        for r in new_rows[1:]:
            if float(np.linalg.norm(r - slim[-1])) < min_seg:
                slim[-1] = r.copy()
            else:
                slim.append(r.copy())

        new_pts = np.stack(slim, axis=0)
        if new_pts.shape[0] < 2 or open_polyline_min_segment_length(new_pts) < _MIN_SEGMENT_LENGTH_M:
            _queue_web_nodes(pts)
            out.append(st)
            continue
        _queue_web_nodes(new_pts)
        out.append(replace(st, midline_twisted_pitch_axis_blade_frame=new_pts))

    # Web half of each junction, applied only now: detection above ran against
    # the original webs, and an inserted point already lies on the edge it is
    # inserted into, so no web moves and no crossing is invalidated.
    for w_idx, insertions in web_insertions.items():
        web_st = out[w_idx]
        w_pts = np.asarray(web_st.midline_twisted_pitch_axis_blade_frame, dtype=float)
        new_w = _insert_points_on_segments(
            w_pts,
            insertions,
            vertex_coincide_tol_m=vertex_coincide_tol_m,
            min_seg=min_seg,
        )
        if new_w is not None:
            out[w_idx] = replace(web_st, midline_twisted_pitch_axis_blade_frame=new_w)

    return tuple(out)


__all__ = [
    "_CAP_WEB_JUNCTION_MERGE_ARC_FRAC",
    "insert_cap_web_interior_junctions",
]
