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


TE_OPEN_BAND = 0.005
TE_OPEN_OVER_C = 0.005
TE_OPEN_OVER_T = 2.0


def _te_frac_drop(cut) -> float:
    """Deprecated: t/c inset that invented wrap on a sharp OML. Prefer ``te_is_blunt``."""
    tc = _cut_tc(cut)
    if tc is None:
        tc = 0.2
    return float(TE_ARC) * float(np.clip(float(tc) / 0.45, 0.25, 1.0))


def te_blunt_eps(cut) -> float:
    """Opening above ``2 t + 0.005 c`` is a finite TE edge, not a point TE plus skin offset."""
    c = max(float(getattr(cut, "chord", 0.0) or 0.0), 1e-12)
    ts = [
        float(w.thickness)
        for w in getattr(cut, "walls", ())
        if getattr(w, "kind", "") in ("skin", "web") and float(getattr(w, "thickness", 0.0) or 0.0) > 0.0
    ]
    t = min(ts) if ts else 0.0
    return TE_OPEN_OVER_T * t + TE_OPEN_OVER_C * c


def te_opening(cut, pts: np.ndarray) -> float:
    """Distance between suction and pressure TE corners on the skin-mid loop."""
    if pts is None or len(pts) < 2:
        return 0.0
    _i_le, i_su, i_pr = _te_corner_indices(pts, cut)
    return float(np.linalg.norm(np.asarray(pts[i_su], float)[:2] - np.asarray(pts[i_pr], float)[:2]))


def te_is_blunt(cut, pts: np.ndarray) -> bool:
    """True when the midline has a finite TE edge, not a single ridge vertex."""
    return te_opening(cut, pts) > te_blunt_eps(cut)


def te_blunt_of_cut(cut) -> bool:
    """``te_is_blunt`` on the cut's skin-mid loop."""
    z = float(getattr(cut, "z", 0.0))
    return te_is_blunt(cut, _skin_mid_loop(cut.shell_graph(), z))


def hold_te_blunt(openings, eps) -> list[bool]:
    """Hysteresis along a span sample: leave blunt when opening ≤ ε; enter when > ε."""
    openings = [float(x) for x in openings]
    eps = [float(x) for x in np.atleast_1d(eps)]
    if len(eps) == 1 and len(openings) > 1:
        eps = eps * len(openings)
    if not openings:
        return []
    state = openings[0] > eps[0]
    out: list[bool] = []
    for o, e in zip(openings, eps):
        e = max(float(e), 1e-15)
        if state and o <= e:
            state = False
        elif (not state) and o > e:
            state = True
        out.append(state)
    return out


def te_blunt_flags_for_cuts(cuts) -> list[bool | None]:
    """Held blunt flags for a station sequence; ``None`` cuts stay ``None``."""
    openings: list[float] = []
    eps: list[float] = []
    idx: list[int] = []
    for i, cut in enumerate(cuts):
        if cut is None:
            continue
        pts = _skin_mid_loop(cut.shell_graph(), float(cut.z))
        openings.append(te_opening(cut, pts))
        eps.append(te_blunt_eps(cut))
        idx.append(i)
    held = hold_te_blunt(openings, eps)
    out: list[bool | None] = [None] * len(cuts)
    for i, flag in zip(idx, held):
        out[i] = flag
    return out


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
    """LE index and suction / pressure TE corners on the skin-mid loop.

    Corners are the farthest pair among vertices within ``TE_OPEN_BAND`` of the
    trailing extreme (same order as ``TE_OPEN_OVER_C``). A sharp NACA collapses
    to one ridge; a cylinder / flatback keeps two.
    """
    frac = _chord_fracs(cut, pts)
    i_le = int(np.argmin(frac))
    fmax = float(frac.max())
    idx = np.flatnonzero(frac >= fmax - TE_OPEN_BAND)
    if len(idx) < 2:
        i_te = int(np.argmax(frac))
        return i_le, i_te, i_te
    xy = np.asarray(pts, float)[idx, :2]
    delta = xy[:, None, :] - xy[None, :, :]
    d2 = np.sum(delta * delta, axis=2)
    k = int(np.argmax(d2))
    a, b = divmod(k, len(idx))
    i0, i1 = int(idx[a]), int(idx[b])
    i_ridge = int(np.argmax(frac))
    n = len(pts)

    def _side_idx(step: int) -> set[int]:
        out: set[int] = set()
        i = i_ridge
        for _ in range(n):
            i = (i + int(step)) % n
            if i == i_le:
                break
            out.add(i)
        return out

    cw, ccw = _side_idx(1), _side_idx(-1)
    same = (i0 in cw and i1 in cw) or (i0 in ccw and i1 in ccw)
    if same:
        return i_le, i_ridge, i_ridge
    if float(pts[i0, 1]) >= float(pts[i1, 1]):
        return i_le, i0, i1
    return i_le, i1, i0


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


def _side_arcs(pts: np.ndarray, cut, i_le: int | None = None, i_su: int | None = None,
               i_pr: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Suction and pressure arcs, each TE-corner → LE (open; TE corners may differ)."""
    if i_le is None or i_su is None or i_pr is None:
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


