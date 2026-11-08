"""Station walks and TE wrap for mapped loft."""
from __future__ import annotations

import numpy as np

from implicit_geometry_kernel.blade.shell_graph import ShellGraph
from implicit_geometry_kernel.blade.topology.graph import (
    TE_ARC,
    _cut_tc,
)


def _chord_fracs(cut, xy: np.ndarray) -> np.ndarray:
    t = np.radians(float(cut.twist))
    c, s = np.cos(t), np.sin(t)
    x = xy[:, 0] - float(cut.offset[0])
    y = xy[:, 1] - float(cut.offset[1])
    local_x = c * x - s * y
    return (local_x + float(cut.pitch_axis) * float(cut.chord)) / max(float(cut.chord), 1e-12)


def _as_xyz(pts: np.ndarray, z: float) -> np.ndarray:
    pts = np.asarray(pts, float)
    if pts.ndim != 2 or len(pts) == 0:
        return np.zeros((0, 3))
    if pts.shape[1] == 2:
        return np.column_stack([pts, np.full(len(pts), float(z))])
    return pts[:, :3]


def _skin_mid_loop(graph: ShellGraph, z: float) -> np.ndarray:
    """Closed skin-mid polyline chained by shared endpoints (closing vertex omitted)."""
    chain = graph.closed_loop_chain(rel_tol=1e-9) if hasattr(graph, "closed_loop_chain") else None
    if chain:
        by_name = {r.name: r for r in graph.loop}
        chunks: list[np.ndarray] = []
        for name, _xy, rev in chain:
            run = by_name.get(name)
            if run is None:
                continue
            pts = _as_xyz(np.asarray(run.points, float), z)
            if rev:
                pts = pts[::-1]
            chunks.append(pts)
        if chunks:
            pts = chunks[0]
            for nxt in chunks[1:]:
                if np.allclose(pts[-1], nxt[0]):
                    pts = np.vstack([pts, nxt[1:]])
                else:
                    pts = np.vstack([pts, nxt])
            if len(pts) > 2 and np.allclose(pts[0], pts[-1]):
                pts = pts[:-1]
            return pts
    try:
        xy = np.asarray(graph.loop_polygon(), float)
    except (ValueError, AttributeError):
        xy = np.zeros((0, 2))
    if len(xy) == 0:
        return np.zeros((0, 3))
    if len(xy) > 2 and np.allclose(xy[0], xy[-1]):
        xy = xy[:-1]
    return np.column_stack([xy[:, :2], np.full(len(xy), float(z))])


def _te_frac_drop(cut) -> float:
    """Chord-fraction length of the TE wrap: full ``TE_ARC`` on a cylinder, a
    fraction of that on a thin foil, so the cap width follows t/c instead of
    switching on a y-span threshold (that switch wrinkled the loft)."""
    tc = _cut_tc(cut)
    if tc is None:
        tc = 0.2
    return float(TE_ARC) * float(np.clip(float(tc) / 0.45, 0.25, 1.0))


def _oml_te_xyz(cut) -> np.ndarray | None:
    """Trailing-edge point of the skin walls (nearest stand-in for OML TE)."""
    key_pts = getattr(cut, "key_points", None)
    if callable(key_pts):
        te = key_pts().get("TE")
        if te is not None:
            return np.asarray(te, float)[:3]
    skins = getattr(cut, "walls_of", lambda *_: ())("skin")
    if not skins:
        return None
    oml = np.vstack([np.asarray(w.points, float) for w in skins])
    fr = _chord_fracs(cut, oml)
    return oml[int(np.argmax(fr))].copy()


def _te_corner_indices(pts: np.ndarray, cut) -> tuple[int, int, int]:
    """LE and the OML-TE ridge index (both corners). Side arcs start at the ridge;
    wrap corners are interpolated along those arcs at ``_te_frac_drop``.
    """
    frac = _chord_fracs(cut, pts)
    i_le = int(np.argmin(frac))
    oml = _oml_te_xyz(cut)
    if oml is not None:
        i_te = int(np.argmin(np.linalg.norm(pts[:, :2] - np.asarray(oml, float)[:2], axis=1)))
        if i_te != i_le:
            return i_le, i_te, i_te
    i_te = int(np.argmax(frac))
    return i_le, i_te, i_te


def _te_wrap_arc(pts: np.ndarray, a: np.ndarray, b: np.ndarray, i_le: int,
                 i_te: int | None = None) -> np.ndarray | None:
    """Skin-mid walk from TE suction to TE pressure that does not pass the LE.

    The walk must include the OML TE ridge when ``i_te`` is given so the loft
    keys sit on the trailing edge rather than a 0.65-chord shortcut.
    """
    if len(pts) < 2:
        return None
    ia = int(np.argmin(np.linalg.norm(pts - np.asarray(a, float), axis=1)))
    ib = int(np.argmin(np.linalg.norm(pts - np.asarray(b, float), axis=1)))
    le = pts[int(i_le)]
    ridge_i = None if i_te is None else int(i_te)

    def has_le(arc: np.ndarray) -> bool:
        return any(np.allclose(p, le) for p in arc)

    def has_ridge(arc: np.ndarray) -> bool:
        if ridge_i is None or ridge_i < 0 or ridge_i >= len(pts):
            return True
        ridge = pts[ridge_i]
        return any(np.allclose(p, ridge) for p in arc)

    def length(arc: np.ndarray) -> float:
        return float(np.sum(np.linalg.norm(np.diff(arc, axis=0), axis=1)))

    if ia == ib:
        if ridge_i is None or ridge_i == ia:
            return None
        via = np.vstack((pts[ia], pts[ridge_i], pts[ib]))
        return via if not has_le(via) else None

    cands = [_walk_arc_dir(pts, ia, ib, 1), _walk_arc_dir(pts, ia, ib, -1)]
    with_ridge = [arc for arc in cands if len(arc) > 1 and not has_le(arc) and has_ridge(arc)]
    if with_ridge:
        return min(with_ridge, key=length)
    if ridge_i is not None and ridge_i != ia and ridge_i != ib:
        via = np.vstack((_walk_arc(pts, ia, ridge_i), _walk_arc(pts, ridge_i, ib)[1:]))
        if len(via) > 1 and not has_le(via):
            return via
    ok = [arc for arc in cands if len(arc) > 1 and not has_le(arc)]
    if not ok:
        return None
    return min(ok, key=length)


def _orient_loop_te(pts: np.ndarray, cut) -> np.ndarray:
    if len(pts) == 0:
        return pts
    _i_le, i_su, _i_pr = _te_corner_indices(pts, cut)
    mid = 0.5 * (pts[i_su] + pts[_i_pr])
    i = int(np.argmin(np.linalg.norm(pts - mid, axis=1)))
    return np.roll(pts, -i, axis=0)


def _walk_arc_dir(pts: np.ndarray, i0: int, i1: int, step: int) -> np.ndarray:
    n = len(pts)
    if n == 0:
        return pts
    out = [pts[i0]]
    if i0 == i1:
        return np.asarray(out, float)
    i = i0
    for _ in range(n):
        i = (i + int(step)) % n
        out.append(pts[i])
        if i == i1:
            break
    return np.asarray(out, float)


def _walk_arc(pts: np.ndarray, i0: int, i1: int) -> np.ndarray:
    n = len(pts)
    if n == 0:
        return pts
    if i0 == i1:
        return pts[i0 : i0 + 1]
    d1, d2 = (i1 - i0) % n, (i0 - i1) % n
    step = 1 if d1 <= d2 else -1
    return _walk_arc_dir(pts, i0, i1, step)


def _xyz_at_frac(arc: np.ndarray, fracs: np.ndarray, f: float) -> np.ndarray:
    """Point on an open arc whose chord fraction is ``f`` (linear between vertices)."""
    if len(arc) == 0:
        return np.zeros(3)
    f = float(np.clip(f, 0.0, 1.0))
    fracs = np.asarray(fracs, float)
    for i in range(len(arc) - 1):
        a, b = float(fracs[i]), float(fracs[i + 1])
        if (a - f) * (b - f) <= 0.0 and abs(b - a) > 1e-15:
            u = (f - a) / (b - a)
            return arc[i] + u * (arc[i + 1] - arc[i])
    return arc[int(np.argmin(np.abs(fracs - f)))].copy()


def _side_arcs(pts: np.ndarray, cut) -> tuple[np.ndarray, np.ndarray]:
    """Suction and pressure arcs, each TE-corner → LE (open; TE corners may differ)."""
    i_le, i_su, i_pr = _te_corner_indices(pts, cut)
    n = len(pts)
    if i_su == i_pr:
        cw = _walk_arc_dir(pts, i_su, i_le, 1)
        ccw = _walk_arc_dir(pts, i_su, i_le, -1)
        y0 = float(cut.offset[1])
        interior = cw[1:-1] if len(cw) > 2 else cw
        if float(np.mean(interior[:, 1])) >= y0:
            return cw, ccw
        return ccw, cw

    def _away(i_from: int, i_other: int) -> int:
        if (i_from + 1) % n == i_other:
            return -1
        if (i_from - 1) % n == i_other:
            return 1
        d_le = (i_le - i_from) % n
        return 1 if d_le <= n - d_le else -1

    return (
        _walk_arc_dir(pts, i_su, i_le, _away(i_su, i_pr)),
        _walk_arc_dir(pts, i_pr, i_le, _away(i_pr, i_su)),
    )


