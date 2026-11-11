"""Rib faces from mapped loft samples."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from implicit_geometry_kernel.blade.stations import TransverseMember
from implicit_geometry_kernel.blade.topology.graph import WallRegionGraph, _cut_tc
from implicit_geometry_kernel.blade.topology.map import (
    _key_name,
    sample_region_polylines,
)
from implicit_geometry_kernel.blade.topology.map_loft import (
    _append_ring,
    _web_pts_from_su,
)
from implicit_geometry_kernel.blade.topology.map_walk import (
    _chord_fracs,
    _skin_mid_loop,
    te_blunt_flags_for_cuts,
    te_is_blunt,
)
from implicit_geometry_kernel.blade.topology.template import ShellTopologyTemplate

def _web_cell_end(vid, web_keys: set[str]) -> bool:
    name = _key_name(vid)
    return name in web_keys or name == "te"


def _group_skin_branches(branches, web_keys: set[str]) -> list[list]:
    """LE→TE groups split only at live web feet (cap anchors stay on the skin)."""
    groups: list[list] = []
    cur: list = []
    for branch in branches:
        cur.append(branch)
        if _web_cell_end(branch.end, web_keys):
            groups.append(cur)
            cur = []
    if cur:
        groups.append(cur)
    return groups


def _concat_branch_pts(branches, samples: Mapping[str, np.ndarray]) -> np.ndarray | None:
    parts: list[np.ndarray] = []
    for branch in branches:
        if branch.name not in samples:
            continue
        pts = np.asarray(samples[branch.name], float)
        if len(pts) < 2:
            continue
        _append_ring(parts, pts)
    if not parts:
        return None
    return np.vstack(parts)


def _close_rib_ring(ring: np.ndarray) -> np.ndarray | None:
    kept = [ring[0]]
    for p in ring[1:]:
        if not np.allclose(kept[-1], p, atol=1e-12, rtol=0.0):
            kept.append(p)
    ring = np.asarray(kept, float)
    if len(ring) > 2 and np.allclose(ring[0], ring[-1], atol=1e-12, rtol=0.0):
        ring = ring[:-1]
    if len(ring) < 3:
        return None
    return ring


def rib_cell_polygons_from_samples(
    template: ShellTopologyTemplate,
    samples: Mapping[str, np.ndarray],
    member: TransverseMember | str,
) -> list[tuple[str, np.ndarray]]:
    """One n-gon per web-bounded cell from loft skin/web nodes (no Steiner vertices).

    Cap ``LoopBranch`` breakpoints stay on the skin midsurface; they do not
    split rib patches. Patch count is ``n_live_webs + 1``.
    """
    label = member if isinstance(member, str) else str(member.name)
    su = [b for b in template.loop if b.side == "suction" and not b.wrap]
    pr = [b for b in template.loop if b.side == "pressure" and not b.wrap]
    wrap = next((b for b in template.loop if b.wrap), None)
    webs = list(template.webs)
    by_key = {_key_name(w.suction): w for w in webs}
    web_keys = set(by_key)
    su_groups = _group_skin_branches(su, web_keys)
    pr_groups = _group_skin_branches(pr, web_keys)
    n_cells = min(len(su_groups), len(pr_groups))
    out: list[tuple[str, np.ndarray]] = []
    for k, (su_g, pr_g) in enumerate(zip(su_groups, pr_groups)):
        su_pts = _concat_branch_pts(su_g, samples)
        pr_pts = _concat_branch_pts(pr_g, samples)
        if su_pts is None or pr_pts is None or len(su_pts) < 2 or len(pr_pts) < 2:
            continue
        aft = by_key.get(_key_name(su_g[-1].end))
        fwd = by_key.get(_key_name(su_g[0].start))
        aft_pts = _web_pts_from_su(samples, aft, su_pts[-1])
        fwd_pts = _web_pts_from_su(samples, fwd, su_pts[0])
        wrap_pts = None
        if wrap is not None and wrap.name in samples and _key_name(su_g[-1].end) == "te":
            wrap_pts = np.asarray(samples[wrap.name], float)
            if len(wrap_pts) and np.linalg.norm(wrap_pts[0] - su_pts[-1]) > np.linalg.norm(
                wrap_pts[-1] - su_pts[-1]
            ):
                wrap_pts = wrap_pts[::-1].copy()
        parts: list[np.ndarray] = []
        _append_ring(parts, su_pts)
        _append_ring(parts, aft_pts)
        _append_ring(parts, wrap_pts)
        _append_ring(parts, pr_pts[::-1])
        if fwd_pts is not None:
            _append_ring(parts, fwd_pts[::-1])
        if not parts:
            continue
        ring = _close_rib_ring(np.vstack(parts))
        if ring is None:
            continue
        name = label if n_cells == 1 else f"{label}:cell_{k}"
        out.append((name, ring))
    return out


def _ring_basis(ring: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Best-fit plane of a 3-D n-gon: origin, u, v, and 2-D xy of the ring."""
    pts = np.asarray(ring, float)
    origin = pts.mean(axis=0)
    rel = pts - origin
    n = np.zeros(3, float)
    for i in range(len(pts)):
        a, b = pts[i], pts[(i + 1) % len(pts)]
        n[0] += (a[1] - b[1]) * (a[2] + b[2])
        n[1] += (a[2] - b[2]) * (a[0] + b[0])
        n[2] += (a[0] - b[0]) * (a[1] + b[1])
    ln = float(np.linalg.norm(n))
    if ln < 1e-15:
        _, _, vh = np.linalg.svd(rel, full_matrices=False)
        n = vh[-1]
        ln = float(np.linalg.norm(n)) or 1.0
    n = n / ln
    seed = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(n, seed)
    u = u / (float(np.linalg.norm(u)) or 1.0)
    v = np.cross(n, u)
    xy = np.column_stack([rel @ u, rel @ v])
    return origin, u, v, xy


def _project_xy(pts: np.ndarray, origin: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    rel = np.asarray(pts, float) - origin
    return np.column_stack([rel @ u, rel @ v])


def _area2(xy: np.ndarray) -> float:
    x, y = xy[:, 0], xy[:, 1]
    return float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _cross2(a, b, c) -> float:
    return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def _in_tri(p, a, b, c, *, eps: float = 1e-12) -> bool:
    c1, c2, c3 = _cross2(a, b, p), _cross2(b, c, p), _cross2(c, a, p)
    return not ((c1 < -eps or c2 < -eps or c3 < -eps) and (c1 > eps or c2 > eps or c3 > eps))


def _point_in_ring(p, xy: np.ndarray) -> bool:
    """Even-odd test in the ring plane (includes boundary)."""
    x, y = float(p[0]), float(p[1])
    inside = False
    n = len(xy)
    for i in range(n):
        x1, y1 = xy[i]
        x2, y2 = xy[(i + 1) % n]
        if abs(_cross2(xy[i], xy[(i + 1) % n], p)) <= 1e-12:
            if min(x1, x2) - 1e-12 <= x <= max(x1, x2) + 1e-12 and (
                min(y1, y2) - 1e-12 <= y <= max(y1, y2) + 1e-12
            ):
                return True
        cond = ((y1 > y) != (y2 > y)) and (
            x < (x2 - x1) * (y - y1) / (y2 - y1 + 1e-30) + x1
        )
        if cond:
            inside = not inside
    return inside


def ear_clip_ring(ring: np.ndarray) -> list[np.ndarray]:
    """Triangles that fill ``ring`` using only its vertices (no Steiner points)."""
    pts = np.asarray(ring, float)
    if len(pts) > 2 and np.allclose(pts[0], pts[-1], atol=1e-12, rtol=0.0):
        pts = pts[:-1]
    if len(pts) < 3:
        return []
    if len(pts) == 3:
        return [pts.copy()]
    _origin, _u, _v, xy = _ring_basis(pts)
    if _area2(xy) < 0.0:
        pts = pts[::-1].copy()
        xy = xy[::-1].copy()
    idx = list(range(len(pts)))
    faces: list[np.ndarray] = []
    guard = 0
    while len(idx) > 3 and guard < 8 * len(pts):
        guard += 1
        n = len(idx)
        clipped = False
        for k in range(n):
            i0, i1, i2 = idx[(k - 1) % n], idx[k], idx[(k + 1) % n]
            a, b, c = xy[i0], xy[i1], xy[i2]
            if _cross2(a, b, c) <= 1e-14:
                continue
            if any(
                _in_tri(xy[j], a, b, c)
                for j in idx
                if j not in (i0, i1, i2)
            ):
                continue
            faces.append(np.stack([pts[i0], pts[i1], pts[i2]]))
            del idx[k]
            clipped = True
            break
        if not clipped:
            break
    if len(idx) == 3:
        faces.append(np.stack([pts[idx[0]], pts[idx[1]], pts[idx[2]]]))
    return faces


def rib_cell_triangles_from_samples(
    template: ShellTopologyTemplate,
    samples: Mapping[str, np.ndarray],
    member: TransverseMember | str,
) -> list[tuple[str, np.ndarray]]:
    """Ear-clip of :func:`rib_cell_polygons_from_samples` (drawing / STL only)."""
    out: list[tuple[str, np.ndarray]] = []
    for name, ring in rib_cell_polygons_from_samples(template, samples, member):
        for tri in ear_clip_ring(ring):
            out.append((name, tri))
    return out


def rib_face_polygons(geo, *, n_chord: int = 12) -> list[tuple[str, np.ndarray]]:
    """Rib n-gons at each ``TransverseMember.z``, vertices from the mapped loft.

    The tip closure is not a rib: see :func:`tip_closure_polygons`.
    """
    out: list[tuple[str, np.ndarray]] = []
    graph = getattr(geo, "region_graph", None)
    if graph is None:
        return out
    for member in getattr(geo, "transverse", ()) or ():
        out.extend(_face_polygons_at(geo, graph, float(member.z), member, int(n_chord)))
    return out


def tip_closure_polygons(geo, *, n_chord: int = 12) -> list[tuple[str, np.ndarray]]:
    """The tip closure's n-gons on the tip station, vertices from the mapped loft.

    Its perimeter vertices are the skin (and any live web) samples of the tip
    station, so it closes the skin mid-surface loop. One n-gon per web-bounded
    cell, named ``tip_closure`` or ``tip_closure:cell_k``.
    """
    graph = getattr(geo, "region_graph", None)
    tip = getattr(geo, "tip_closure", None)
    if graph is None or tip is None:
        return []
    return _face_polygons_at(geo, graph, float(geo.tip_closure_z), str(tip.name), int(n_chord))


def _face_polygons_at(geo, graph, z: float, member, n_chord: int) -> list[tuple[str, np.ndarray]]:
    cut = section_cut_near(geo, z)
    if cut is None:
        return []
    sg = cut.shell_graph()
    tmpl = _template_from_cut(graph, cut, z, n_chord)
    samples = sample_region_polylines(cut, sg, tmpl, graph)
    return rib_cell_polygons_from_samples(tmpl, samples, member)


def rib_face_triangles(geo, *, n_chord: int = 12) -> list[tuple[str, np.ndarray]]:
    """STL tessellation: ear-clip of each rib n-gon (same vertices, no centroid)."""
    faces = []
    for name, ring in rib_face_polygons(geo, n_chord=n_chord):
        for tri in ear_clip_ring(ring):
            faces.append((name, tri))
    return faces


def section_cut_near(geo, z: float):
    """``section_at(z)``, or a neighbour if that station is infeasible."""
    from implicit_geometry_kernel.contract import InfeasibleSectionGeometry

    z = float(z)
    try:
        return geo.section_at(z)
    except InfeasibleSectionGeometry:
        pass
    for dz in (1e-6, -1e-6, 1e-4, -1e-4, 1e-3, -1e-3):
        try:
            return geo.section_at(z + dz)
        except (InfeasibleSectionGeometry, ValueError):
            continue
    return None


def _live_web_fracs(cut) -> list[float]:
    """Twist-aware live web chord fractions from the shell graph."""
    live = []
    for run in cut.shell_graph().webs:
        pts = np.asarray(run.points, float)[:, :2]
        live.append(float(np.clip(_chord_fracs(cut, pts[:1])[0], 0.0, 1.0)))
    return live


def _template_from_cut(graph: WallRegionGraph, cut, z: float, n_chord: int,
                       te_blunt: bool | None = None) -> ShellTopologyTemplate:
    dropped = getattr(cut, "dropped_web_indices", None)
    shell = cut.shell_graph()
    if te_blunt is None:
        pts = _skin_mid_loop(shell, float(z))
        te_blunt = te_is_blunt(cut, pts)
    return ShellTopologyTemplate.from_region_graph(
        graph, float(z), n_chord=int(n_chord),
        live_web_fracs=_live_web_fracs(cut) or None, t_c=_cut_tc(cut),
        dropped_web_indices=None if dropped is None else tuple(int(i) for i in dropped),
        te_blunt=te_blunt,
        side_fracs=_station_side_fracs(graph, cut, shell, float(z), te_blunt),
    )


def _station_side_fracs(graph: WallRegionGraph, cut, shell, z: float,
                        te_blunt: bool) -> dict[str, dict[str, float]]:
    """Each side's chord fraction of every anchor key, where the sampler puts it."""
    from implicit_geometry_kernel.blade.topology.map import place_station_keys

    keys = place_station_keys(cut, shell, graph, z, te_blunt=te_blunt)
    out: dict[str, dict[str, float]] = {"suction": {}, "pressure": {}}
    for (name, side), xyz in keys.items():
        if name in ("le", "te") or side not in out:
            continue
        out[side][str(name)] = float(_chord_fracs(cut, np.asarray(xyz, float)[None, :2])[0])
    return out


def templates_along(geo, zs: Sequence[float], n_chord: int = 12) -> list[ShellTopologyTemplate]:
    graph = geo.region_graph
    cuts = [section_cut_near(geo, float(z)) for z in zs]
    flags = te_blunt_flags_for_cuts(cuts)
    out = []
    for cut, flag in zip(cuts, flags):
        if cut is None:
            continue
        out.append(_template_from_cut(graph, cut, float(cut.z), n_chord, te_blunt=flag))
    return out
