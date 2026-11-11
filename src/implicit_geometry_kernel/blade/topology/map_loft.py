"""Loft mapped strips between sampled stations."""
from __future__ import annotations

from typing import Mapping

import numpy as np

from implicit_geometry_kernel.blade.loft.locked_loft import arc_resample_open
from implicit_geometry_kernel.blade.topology.sizing import resample_open_h
from implicit_geometry_kernel.blade.topology.template import (
    LoopBranch,
    ShellTopologyTemplate,
    WebBranch,
)

def loft_mapped_strips(s0: Mapping[str, np.ndarray], s1: Mapping[str, np.ndarray],
                       t0: ShellTopologyTemplate, t1: ShellTopologyTemplate
                       ) -> list[tuple[str, np.ndarray, bool]]:
    """All-quad strips between two sampled stations. Omit one-sided webs.

    Each t1 skin branch is ruled from the t0 skin walk between its mapped
    endpoints (web onset/drop splits a panel along the mould, not a chord).
    Both sides of every strip are the stations' own nodes, so neighbouring
    bands share every station node.

    Where any strip of the band changes node count (a web or cap starts or
    ends, a panel is born from a point, the TE goes blunt or sharp), the band
    gets a mid-span row: equal strips split into two quads, and the others take
    an all-quad transition (:func:`_transition_quads`).

    Returns ``(name, face, collapsing)``; ``face`` is a quad ``(4, 3)``.
    ``collapsing`` marks the faces of a transition (and the triangle halves of
    a folded quad, which are ``(3, 3)``).
    """
    pairs = _band_pairs(s0, s1, t0, t1)
    event = any(len(pa) != len(pb) for _name, pa, pb in pairs)
    faces: list[tuple[str, np.ndarray, bool]] = []
    for name, pa, pb in pairs:
        if not event:
            faces.extend(_ruled_quads(name, pa, pb))
        elif len(pa) == len(pb):
            faces.extend(_split_quads(name, pa, pb))
        else:
            faces.extend(_transition_quads(name, pa, pb))
    return faces


def _band_pairs(s0: Mapping[str, np.ndarray], s1: Mapping[str, np.ndarray],
                t0: ShellTopologyTemplate, t1: ShellTopologyTemplate
                ) -> list[tuple[str, np.ndarray, np.ndarray]]:
    """``(name, t0 nodes, t1 nodes)`` for every strip of the band."""
    raw: list[tuple[str, np.ndarray, np.ndarray]] = []
    t0_wraps = any(b.wrap for b in t0.loop)
    for branch in t1.loop:
        if branch.wrap and not t0_wraps and branch.name in s1:
            ridge0 = _sample_te_ridge(s0)
            if ridge0 is not None:
                raw.append((branch.name, ridge0.reshape(1, -1), np.asarray(s1[branch.name], float)))
            continue
        if branch.wrap or branch.name not in s1:
            continue
        pb = np.asarray(s1[branch.name], float)
        pa = _skin_parent(branch, s0, t0, t1, pb)
        if pa is not None:
            raw.append((branch.name, pa, pb))
    for branch in t0.loop:
        if not branch.wrap or branch.name not in s0:
            continue
        dest = _strip_dest(branch, s1, t0, t1)
        pa = np.asarray(s0[branch.name], float)
        if dest is not None and dest in s1:
            pb = np.asarray(s1[dest], float)
        else:
            ridge = _sample_te_ridge(s1)
            if ridge is None:
                continue
            pb, dest = ridge.reshape(1, -1), branch.name
        raw.append((dest, pa, pb))
    names1 = {w.name for w in t1.webs}
    for web in t0.webs:
        if web.name not in s0 or web.name not in s1 or web.name not in names1:
            continue
        raw.append((web.name, np.asarray(s0[web.name], float), np.asarray(s1[web.name], float)))
    out: list[tuple[str, np.ndarray, np.ndarray]] = []
    for name, pa, pb in raw:
        if len(pa) == 0 or len(pb) == 0 or not (np.all(np.isfinite(pa)) and np.all(np.isfinite(pb))):
            continue
        pa, pb = _dedupe_open(np.asarray(pa, float)), _dedupe_open(np.asarray(pb, float))
        if len(pa) == 1 and len(pb) == 1:
            continue
        out.append((name, pa, pb))
    return out


def _ruled_quads(name: str, pa: np.ndarray, pb: np.ndarray) -> list[tuple[str, np.ndarray, bool]]:
    out: list[tuple[str, np.ndarray, bool]] = []
    for i in range(len(pa) - 1):
        quad = np.stack([pa[i], pa[i + 1], pb[i + 1], pb[i]])
        if _quad_bowtie(quad):
            a, b, c, d = quad
            _append_face(out, name, np.stack([a, b, c]), True)
            _append_face(out, name, np.stack([a, c, d]), True)
        else:
            _append_face(out, name, quad, False)
    return out


def _mid(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Midpoint of a spanwise or interior edge, ``p`` on t0 first (same bits on both sides)."""
    return 0.5 * (p + q)


def _split_quads(name: str, pa: np.ndarray, pb: np.ndarray) -> list[tuple[str, np.ndarray, bool]]:
    """Two quads per ruled quad, through the mid-span row of an event band."""
    m = [_mid(pa[i], pb[i]) for i in range(len(pa))]
    out: list[tuple[str, np.ndarray, bool]] = []
    for i in range(len(pa) - 1):
        _append_face(out, name, np.stack([pa[i], pa[i + 1], m[i + 1], m[i]]), False)
        _append_face(out, name, np.stack([m[i], m[i + 1], pb[i + 1], pb[i]]), False)
    return out


def _transition_quads(name: str, pa: np.ndarray, pb: np.ndarray) -> list[tuple[str, np.ndarray, bool]]:
    """All-quad strip between station polylines with different node counts.

    Both sides have an even number of edges (each region has an odd node
    count, and onset keys land on even nodes), so every other node is a
    coarse node. The coarse sides are zipped into triangles and each triangle
    splits into three quads through its centroid and edge midpoints. On a
    station line the midpoint of a coarse edge is the station node between
    its ends, and on the band's spanwise sides it is the mid-span row node,
    so the strip shares every node with its neighbours. A side of one node
    (a panel born from a point) is the apex of a fan.
    """
    if (len(pa) - 1) % 2 or (len(pb) - 1) % 2:
        raise ValueError(
            f"transition strip {name!r} has {len(pa) - 1} and {len(pb) - 1} edges; "
            "an all-quad transition needs an even number on each station"
        )
    ca, cb = pa[::2], pb[::2]

    def point(v: tuple[str, int]) -> np.ndarray:
        return ca[v[1]] if v[0] == "a" else cb[v[1]]

    def mid(u: tuple[str, int], v: tuple[str, int]) -> np.ndarray:
        if u[0] == v[0] and abs(u[1] - v[1]) == 1:
            fine = pa if u[0] == "a" else pb
            return fine[2 * min(u[1], v[1]) + 1]
        if u[0] == "b" and v[0] == "a":
            u, v = v, u
        return _mid(point(u), point(v))

    tris: list[tuple[tuple[str, int], tuple[str, int], tuple[str, int]]] = []
    if len(ca) == 1:
        tris = [(("a", 0), ("b", j + 1), ("b", j)) for j in range(len(cb) - 1)]
    elif len(cb) == 1:
        tris = [(("a", i), ("a", i + 1), ("b", 0)) for i in range(len(ca) - 1)]
    else:
        fa, fb = _arc_fracs(ca), _arc_fracs(cb)
        i = j = 0
        while i < len(ca) - 1 or j < len(cb) - 1:
            if j == len(cb) - 1 or (i < len(ca) - 1 and fa[i + 1] <= fb[j + 1]):
                tris.append((("a", i), ("a", i + 1), ("b", j)))
                i += 1
            else:
                tris.append((("a", i), ("b", j + 1), ("b", j)))
                j += 1
    out: list[tuple[str, np.ndarray, bool]] = []
    for v0, v1, v2 in tris:
        p0, p1, p2 = point(v0), point(v1), point(v2)
        c = (p0 + p1 + p2) / 3.0
        m01, m12, m20 = mid(v0, v1), mid(v1, v2), mid(v2, v0)
        _append_face(out, name, np.stack([p0, m01, c, m20]), True)
        _append_face(out, name, np.stack([p1, m12, c, m01]), True)
        _append_face(out, name, np.stack([p2, m20, c, m12]), True)
    return out


def _dedupe_open(pts: np.ndarray) -> np.ndarray:
    keep = [pts[0]]
    for q in pts[1:]:
        if not np.allclose(keep[-1], q, atol=1e-12, rtol=0.0):
            keep.append(q)
    return np.asarray(keep, float)


def _arc_fracs(pts: np.ndarray) -> np.ndarray:
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    return s / s[-1] if s[-1] > 0.0 else np.linspace(0.0, 1.0, len(pts))


def _match_open_counts(pa: np.ndarray, pb: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Arc-resample both open polylines to the same vertex count (keep endpoints)."""
    pa, pb = np.asarray(pa, float), np.asarray(pb, float)
    n = max(len(pa), len(pb), 2)
    if len(pa) != n:
        pa = arc_resample_open(pa, n)
    if len(pb) != n:
        pb = arc_resample_open(pb, n)
    return pa, pb


def _finite_loft_face(face: np.ndarray) -> bool:
    """Drop collapsed or non-finite strips (they poison 3-D autoscale)."""
    p = np.asarray(face, float)
    if p.ndim != 2 or len(p) < 3 or p.shape[1] < 3:
        return False
    if not np.all(np.isfinite(p)):
        return False
    uniq = [p[0]]
    for q in p[1:]:
        if not np.allclose(uniq[-1], q, atol=1e-12, rtol=0.0):
            uniq.append(q)
    if len(uniq) >= 2 and np.allclose(uniq[0], uniq[-1], atol=1e-12, rtol=0.0):
        uniq = uniq[:-1]
    return len(uniq) >= 3


def _append_face(faces: list[tuple[str, np.ndarray, bool]], name: str,
                 face: np.ndarray, collapsing: bool) -> None:
    if _finite_loft_face(face):
        faces.append((name, np.asarray(face, float), collapsing))


def _xyz_of_vid(samples: Mapping[str, np.ndarray], tmpl: ShellTopologyTemplate, vid) -> np.ndarray | None:
    for b in tmpl.loop:
        if b.name not in samples:
            continue
        pts = np.asarray(samples[b.name], float)
        if len(pts) == 0:
            continue
        if b.start == vid:
            return pts[0].copy()
        if b.end == vid:
            return pts[-1].copy()
    return None


def _concat_side_samples(samples: Mapping[str, np.ndarray], tmpl: ShellTopologyTemplate,
                         side: str) -> np.ndarray | None:
    """LE→TE skin polyline on one side from sampled loop branches."""
    parts: list[np.ndarray] = []
    for b in tmpl.loop:
        if b.wrap or b.side != side or b.name not in samples:
            continue
        _append_ring(parts, np.asarray(samples[b.name], float))
    if not parts:
        return None
    return np.vstack(parts)


def _resample_keys(p0: np.ndarray, p1: np.ndarray, arc: np.ndarray, n: int,
                   *, bunch_end: bool = False, size_field: bool = False,
                   chord: float = 1.0, t: float = 0.0) -> np.ndarray:
    """Arc-length resample between two keys that already lie on ``arc`` (TE→LE)."""
    n = max(int(n), 2)
    if len(arc) < 2:
        out = np.repeat(np.asarray(p0, float).reshape(1, -1), n, axis=0)
        out[-1] = p1
        return out
    i0 = int(np.argmin(np.linalg.norm(arc - p0, axis=1)))
    i1 = int(np.argmin(np.linalg.norm(arc - p1, axis=1)))
    if i0 == i1:
        sub = np.vstack([p0, p1])
    elif i0 < i1:
        sub = np.vstack([p0, arc[i0 + 1 : i1], p1])
    else:
        sub = np.vstack([p0, arc[i0 - 1 : i1 : -1], p1])
    if size_field:
        sampled = resample_open_h(sub, chord, t)
    elif bunch_end and n > 2:
        sampled = _arc_resample_bunched(sub, n)
    else:
        sampled = arc_resample_open(sub, n)
    sampled[0] = np.asarray(p0, float)
    sampled[-1] = np.asarray(p1, float)
    return sampled


def _arc_resample_bunched(pts: np.ndarray, n: int) -> np.ndarray:
    """Arc-length resample clustered toward the last vertex (TE)."""
    pts = np.asarray(pts, float)
    n = max(int(n), 2)
    if len(pts) < 2:
        return arc_resample_open(pts, n)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 1e-15:
        return arc_resample_open(pts, n)
    u = np.sin(np.linspace(0.0, 0.5 * np.pi, n)) ** 2
    t = u * s[-1]
    out = np.empty((n, pts.shape[1]))
    for j, tj in enumerate(t):
        k = int(np.searchsorted(s, tj, side="right") - 1)
        k = min(max(k, 0), len(seg) - 1)
        span = seg[k] if seg[k] > 1e-15 else 1.0
        w = (tj - s[k]) / span
        out[j] = pts[k] + w * (pts[k + 1] - pts[k])
    return out


def _on_chain(chain: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Nearest even node of the t0 side chain to a key that t0 lacks (XY).

    Region boundaries sit on even chain nodes (every region has an even number
    of edges), so a key placed on an even node keeps both strips beside it at
    an even edge count, which the all-quad transition needs.
    """
    pts = np.asarray(chain, float)[::2]
    d = np.linalg.norm(pts[:, :2] - np.asarray(p, float)[:2], axis=1)
    return pts[int(np.argmin(d))].copy()


def _skin_parent(branch: LoopBranch, s0: Mapping[str, np.ndarray],
                 t0: ShellTopologyTemplate, t1: ShellTopologyTemplate, pb: np.ndarray
                 ) -> np.ndarray | None:
    """t0 skin arc between the mapped endpoints of a t1 branch (not a chord).

    A web that only exists on t1 is placed on the t0 mould by nearest XY, so the
    TE panel splits at the new foot instead of collapsing to the TE point.
    """
    pb = np.asarray(pb, float)
    n = max(len(pb), 2)
    chain = _concat_side_samples(s0, t0, branch.side)
    p0 = _te_key_xyz(s0, t0, branch.start, branch.side)
    p1 = _te_key_xyz(s0, t0, branch.end, branch.side)
    if chain is None:
        return _collapsed_parent(branch, s0, t0, t1, n)
    if p0 is None:
        p0 = _on_chain(chain, pb[0])
    if p1 is None:
        p1 = _on_chain(chain, pb[-1])
    own = _chain_nodes_between(chain, p0, p1)
    if own is not None:
        return own
    return _resample_keys(p0, p1, chain, n)


def _chain_nodes_between(chain: np.ndarray, p0: np.ndarray, p1: np.ndarray) -> np.ndarray | None:
    """The t0 station's own nodes from ``p0`` to ``p1``.

    A strip's lower edge must be the t0 station's nodes, or the strip below and
    the strip above meet at different nodes along the station line. Where the
    two stations carry different counts (a topology event), the band takes an
    all-quad transition instead of resampling either side.
    """
    chain = np.asarray(chain, float)
    d0 = np.linalg.norm(chain - np.asarray(p0, float), axis=1)
    d1 = np.linalg.norm(chain - np.asarray(p1, float), axis=1)
    i0, i1 = int(np.argmin(d0)), int(np.argmin(d1))
    if d0[i0] > 1e-12 or d1[i1] > 1e-12:
        return None
    step = 1 if i1 >= i0 else -1
    return chain[i0 : i1 + step if i1 + step >= 0 else None : step].copy()


def _te_key_xyz(samples: Mapping[str, np.ndarray], tmpl: ShellTopologyTemplate, vid, side: str) -> np.ndarray | None:
    """A t1 key on the t0 station, with the TE ridge and TE corner as one key.

    Where the trailing edge goes from blunt to sharp, a sharp branch ends at the
    ridge ``("TE",)`` and the blunt station has a corner ``("TE", side)`` per
    side instead (and the other way round). The strip must end at that corner,
    or it stops one node short of it and leaves a sliver beside the TE strip.
    """
    if vid in tmpl.vertices:
        got = _xyz_of_vid(samples, tmpl, vid)
        if got is not None:
            return got
    if vid == ("TE",):
        return _xyz_of_vid(samples, tmpl, ("TE", side))
    if isinstance(vid, tuple) and len(vid) == 2 and vid[0] == "TE":
        return _xyz_of_vid(samples, tmpl, ("TE",))
    return None


def _collapsed_parent(branch: LoopBranch, s0: Mapping[str, np.ndarray],
                      t0: ShellTopologyTemplate, t1: ShellTopologyTemplate, n: int) -> np.ndarray | None:
    try:
        start0 = t1.map_vertex(branch.start, t0)
        end0 = t1.map_vertex(branch.end, t0)
    except KeyError:
        return None
    p0 = _xyz_of_vid(s0, t0, start0)
    p1 = _xyz_of_vid(s0, t0, end0)
    if p0 is None:
        p0 = p1
    if p1 is None:
        p1 = p0
    if p0 is None or p1 is None:
        return None
    n = max(int(n), 2)
    if np.allclose(p0, p1):
        return np.repeat(np.asarray(p0, float).reshape(1, -1), n, axis=0)
    return np.linspace(p0, p1, n)


def _same_web_set(t0: ShellTopologyTemplate, t1: ShellTopologyTemplate) -> bool:
    s0, s1 = t0.signature, t1.signature
    if isinstance(s0, tuple) and isinstance(s1, tuple) and len(s0) == 2 and len(s1) == 2:
        return s0[1] == s1[1]
    return s0 == s1


def _sample_te_ridge(samples: Mapping[str, np.ndarray]) -> np.ndarray | None:
    for name in ("te_su", "te_pr", "te"):
        pts = np.asarray(samples.get(name, []), float)
        if len(pts):
            return np.asarray(pts[-1], float)
    return None


def _strip_dest(branch: LoopBranch, s1: Mapping[str, np.ndarray],
                t0: ShellTopologyTemplate, t1: ShellTopologyTemplate) -> str | None:
    if _same_web_set(t0, t1) and branch.name in s1:
        return branch.name
    try:
        start_m = t0.map_vertex(branch.start, t1)
        end_m = t0.map_vertex(branch.end, t1)
    except KeyError:
        return branch.name if branch.name in s1 else None
    if start_m == end_m and not branch.wrap:
        return None
    for other in t1.loop:
        if other.name not in s1 or other.side != branch.side:
            continue
        if other.start == start_m and other.end == end_m:
            return other.name
    return None


def _quad_bowtie(quad: np.ndarray) -> bool:
    """True if the two triangles of a loft quad have opposing normals (a fold)."""
    a, b, c, d = np.asarray(quad, float)
    n1 = np.cross(b - a, c - a)
    n2 = np.cross(c - a, d - a)
    if np.linalg.norm(n1) < 1e-18 or np.linalg.norm(n2) < 1e-18:
        return False
    return float(np.dot(n1, n2)) < -0.2 * float(np.linalg.norm(n1) * np.linalg.norm(n2))


def _append_ring(parts: list[np.ndarray], chain: np.ndarray | None) -> None:
    if chain is None:
        return
    c = np.asarray(chain, float)
    if len(c) == 0:
        return
    if parts and np.allclose(parts[-1][-1], c[0], atol=1e-12, rtol=0.0):
        c = c[1:]
    if len(c):
        parts.append(c)


def _web_pts_from_su(samples: Mapping[str, np.ndarray], web: WebBranch | None,
                     su_end: np.ndarray) -> np.ndarray | None:
    if web is None or web.name not in samples:
        return None
    pts = np.asarray(samples[web.name], float)
    if len(pts) < 2:
        return None
    if np.linalg.norm(pts[-1] - su_end) < np.linalg.norm(pts[0] - su_end):
        pts = pts[::-1].copy()
    return pts

