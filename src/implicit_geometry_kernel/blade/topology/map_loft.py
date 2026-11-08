"""Loft mapped strips between sampled stations."""
from __future__ import annotations

from typing import Mapping

import numpy as np

from implicit_geometry_kernel.blade.locked_loft import arc_resample_open
from implicit_geometry_kernel.blade.topology.template import (
    LoopBranch,
    ShellTopologyTemplate,
    WebBranch,
)

def loft_mapped_strips(s0: Mapping[str, np.ndarray], s1: Mapping[str, np.ndarray],
                       t0: ShellTopologyTemplate, t1: ShellTopologyTemplate
                       ) -> list[tuple[str, np.ndarray, bool]]:
    """Quad (or collapsing) strips between two templates. Omit one-sided webs.

    Each t1 skin branch is ruled from the t0 skin walk between its mapped
    endpoints (web onset/drop splits a panel along the mould, not a chord).

    Returns ``(name, face, collapsing)`` with ``face`` shape ``(4, 3)`` or ``(3, 3)``.
    """
    faces: list[tuple[str, np.ndarray, bool]] = []
    used: set[str] = set()
    for branch in t1.loop:
        if branch.wrap or branch.name not in s1:
            continue
        pb = np.asarray(s1[branch.name], float)
        pa = _skin_parent(branch, s0, t0, t1, pb)
        if pa is None:
            continue
        birth = bool(np.allclose(pa[0], pa[-1]))
        faces.extend(_strip_faces(branch.name, pa, pb, t0, t1, branch, birth=birth))
        used.add(branch.name)
    for branch in t0.loop:
        if not branch.wrap or branch.name not in s0:
            continue
        dest = _strip_dest(branch, s1, t0, t1)
        if dest is None:
            continue
        pa, pb = np.asarray(s0[branch.name], float), np.asarray(s1[dest], float)
        faces.extend(_strip_faces(dest, pa, pb, t0, t1, branch))
        used.add(dest)
    names1 = {w.name for w in t1.webs}
    for web in t0.webs:
        if web.name not in s0 or web.name not in s1 or web.name not in names1:
            continue
        pa, pb = np.asarray(s0[web.name], float), np.asarray(s1[web.name], float)
        n = min(len(pa), len(pb))
        if n < 2:
            continue
        for i in range(n - 1):
            quad = np.stack([pa[i], pa[i + 1], pb[i + 1], pb[i]])
            if _quad_bowtie(quad):
                a, b, c, d = quad
                faces.append((web.name, np.stack([a, b, c]), True))
                faces.append((web.name, np.stack([a, c, d]), True))
                continue
            faces.append((web.name, quad, False))
    return faces


def _strip_faces(name: str, pa: np.ndarray, pb: np.ndarray, t0, t1, branch,
                 *, birth: bool = False) -> list[tuple[str, np.ndarray, bool]]:
    pa, pb = np.asarray(pa, float), np.asarray(pb, float)
    if len(pa) != len(pb) or len(pa) < 2:
        n = min(len(pa), len(pb))
        if n < 2:
            return []
        pa, pb = pa[:n], pb[:n]
    collapsing = birth or np.allclose(pa[0], pa[-1]) or np.allclose(pb[0], pb[-1])
    if not birth:
        try:
            start_m = t0.map_vertex(branch.start, t1) if branch.start in t0.vertices else branch.start
            end_m = t0.map_vertex(branch.end, t1) if branch.end in t0.vertices else branch.end
            collapsing = collapsing or (start_m == end_m and not branch.wrap)
        except KeyError:
            pass
    if branch.wrap and np.allclose(pa[0], pa[-1]) and np.allclose(pb[0], pb[-1]):
        return []
    out = []
    if collapsing:
        for i in range(len(pa) - 1):
            out.append((name, np.stack([pa[i], pa[i + 1], pb[min(i + 1, len(pb) - 1)]]), True))
        return out
    for i in range(len(pa) - 1):
        quad = np.stack([pa[i], pa[i + 1], pb[i + 1], pb[i]])
        if _quad_bowtie(quad):
            a, b, c, d = quad
            out.append((name, np.stack([a, b, c]), True))
            out.append((name, np.stack([a, c, d]), True))
        else:
            out.append((name, quad, False))
    return out


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
                   *, bunch_end: bool = False) -> np.ndarray:
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
    if bunch_end and n > 2:
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
    d = np.linalg.norm(np.asarray(chain, float)[:, :2] - np.asarray(p, float)[:2], axis=1)
    return np.asarray(chain, float)[int(np.argmin(d))].copy()


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
    p0 = _xyz_of_vid(s0, t0, branch.start) if branch.start in t0.vertices else None
    p1 = _xyz_of_vid(s0, t0, branch.end) if branch.end in t0.vertices else None
    if chain is None:
        return _collapsed_parent(branch, s0, t0, t1, n)
    if p0 is None:
        p0 = _on_chain(chain, pb[0])
    if p1 is None:
        p1 = _on_chain(chain, pb[-1])
    return _resample_keys(p0, p1, chain, n)


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

