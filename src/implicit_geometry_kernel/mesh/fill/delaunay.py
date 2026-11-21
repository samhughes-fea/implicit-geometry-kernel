"""Delaunay fill of a planar n-gon: hex Steiner lattice, or an unstructured interior.

The ring vertices are the boundary nodes. Qhull when SciPy is present, otherwise
the numba Bowyer–Watson.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np
from numba import njit

from implicit_geometry_kernel.mesh.fill.polygon import _hex_seeds, _ring_basis, points_in_ring


def _area2(xy: np.ndarray) -> float:
    x, y = xy[:, 0], xy[:, 1]
    return float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _median_edge(xy: np.ndarray) -> float:
    d = np.linalg.norm(np.roll(xy, -1, axis=0) - xy, axis=1)
    d = d[d > 1e-14]
    return float(np.median(d)) if len(d) else 0.0


def _target_size(xy: np.ndarray, size: float | None) -> float:
    area = 0.5 * abs(_area2(xy))
    h = float(size) if size is not None else _median_edge(xy)
    if area > 0.0:
        h = max(h, float(np.sqrt(area / 96.0)))
    return h


def _unique_xy(xy: np.ndarray) -> np.ndarray:
    out = [xy[0]]
    for p in xy[1:]:
        if float(np.linalg.norm(p - out[-1])) > 1e-12:
            out.append(p)
    if len(out) > 2 and float(np.linalg.norm(out[0] - out[-1])) < 1e-12:
        out.pop()
    return np.asarray(out, float)


def _fan(pts: np.ndarray) -> list[np.ndarray]:
    if len(pts) < 3:
        return []
    if len(pts) == 3:
        return [pts.copy()]
    return [np.vstack([pts[0], pts[i], pts[i + 1]]) for i in range(1, len(pts) - 1)]


def _has_scipy() -> bool:
    import importlib.util

    return importlib.util.find_spec("scipy") is not None


@njit(cache=True)
def _in_circum(ax, ay, bx, by, cx, cy, px, py) -> bool:
    orient = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    if orient < 0.0:
        bx, by, cx, cy = cx, cy, bx, by
    adx, ady = ax - px, ay - py
    bdx, bdy = bx - px, by - py
    cdx, cdy = cx - px, cy - py
    det = (
        (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy)
        - (bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady)
        + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady)
    )
    return det > 1e-18


@njit(cache=True)
def _bowyer_watson_numba(points: np.ndarray) -> np.ndarray:
    n = points.shape[0]
    if n < 3:
        return np.empty((0, 3), np.int32)
    lo0 = points[:, 0].min()
    lo1 = points[:, 1].min()
    hi0 = points[:, 0].max()
    hi1 = points[:, 1].max()
    dx = hi0 - lo0
    dy = hi1 - lo1
    d = np.sqrt(dx * dx + dy * dy)
    if d < 1e-15:
        d = 1.0
    mid0 = 0.5 * (lo0 + hi0)
    mid1 = 0.5 * (lo1 + hi1)
    verts = np.empty((n + 3, 2), np.float64)
    verts[:n] = points
    verts[n, 0], verts[n, 1] = mid0, mid1 + 20.0 * d
    verts[n + 1, 0], verts[n + 1, 1] = mid0 - 20.0 * d, mid1 - 20.0 * d
    verts[n + 2, 0], verts[n + 2, 1] = mid0 + 20.0 * d, mid1 - 20.0 * d
    max_t = 8 * n + 16
    tris = np.empty((max_t, 3), np.int32)
    n_tris = 1
    tris[0, 0], tris[0, 1], tris[0, 2] = n, n + 1, n + 2
    bad = np.empty(max_t, np.uint8)
    ea = np.empty(3 * max_t, np.int32)
    eb = np.empty(3 * max_t, np.int32)
    for i in range(n):
        px, py = verts[i, 0], verts[i, 1]
        n_bad = 0
        for t in range(n_tris):
            a, b, c = tris[t, 0], tris[t, 1], tris[t, 2]
            hit = _in_circum(
                verts[a, 0], verts[a, 1], verts[b, 0], verts[b, 1],
                verts[c, 0], verts[c, 1], px, py,
            )
            bad[t] = 1 if hit else 0
            if hit:
                n_bad += 1
        n_e = 0
        for t in range(n_tris):
            if bad[t] == 0:
                continue
            a, b, c = tris[t, 0], tris[t, 1], tris[t, 2]
            for u, v in ((a, b), (b, c), (c, a)):
                if u > v:
                    u, v = v, u
                ea[n_e] = u
                eb[n_e] = v
                n_e += 1
        w = 0
        for t in range(n_tris):
            if bad[t] == 0:
                tris[w, 0], tris[w, 1], tris[w, 2] = tris[t, 0], tris[t, 1], tris[t, 2]
                w += 1
        n_tris = w
        for e in range(n_e):
            cnt = 0
            u, v = ea[e], eb[e]
            for f in range(n_e):
                if ea[f] == u and eb[f] == v:
                    cnt += 1
            if cnt != 1:
                continue
            ax, ay = verts[u, 0], verts[u, 1]
            bx, by = verts[v, 0], verts[v, 1]
            if (bx - ax) * (py - ay) - (by - ay) * (px - ax) > 0.0:
                tris[n_tris, 0], tris[n_tris, 1], tris[n_tris, 2] = u, v, i
            else:
                tris[n_tris, 0], tris[n_tris, 1], tris[n_tris, 2] = v, u, i
            n_tris += 1
    keep = np.empty((n_tris, 3), np.int32)
    k = 0
    for t in range(n_tris):
        a, b, c = tris[t, 0], tris[t, 1], tris[t, 2]
        if a >= n or b >= n or c >= n:
            continue
        keep[k, 0], keep[k, 1], keep[k, 2] = a, b, c
        k += 1
    return keep[:k]


def _simplices(cloud: np.ndarray, engine: str | None) -> np.ndarray:
    if engine is None:
        engine = "scipy" if _has_scipy() else "numba"
    if engine == "scipy":
        from scipy.spatial import Delaunay

        return np.asarray(Delaunay(cloud).simplices, np.int32)
    if engine != "numba":
        raise ValueError(f"unknown Delaunay engine {engine!r}")
    return _bowyer_watson_numba(np.asarray(cloud, np.float64))


def _edge_clearance(xy: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Distance from each point to the nearest edge of ``xy``."""
    pts = np.asarray(pts, float).reshape(-1, 2)
    a = np.asarray(xy, float)
    b = np.roll(a, -1, axis=0)
    ab = b - a
    lab2 = np.maximum(np.sum(ab * ab, axis=1), 1e-30)
    ap = pts[:, None, :] - a[None, :, :]
    t = np.clip(np.sum(ap * ab[None, :, :], axis=2) / lab2, 0.0, 1.0)
    proj = a[None, :, :] + t[:, :, None] * ab[None, :, :]
    return np.linalg.norm(pts[:, None, :] - proj, axis=2).min(axis=1)


def _split_boundary(xy: np.ndarray, h: float) -> np.ndarray:
    """Keep every ring vertex and split edges so each piece is about ``h`` long."""
    if h < 1e-15 or len(xy) < 2:
        return np.asarray(xy, float)
    parts: list[np.ndarray] = []
    n = len(xy)
    for i in range(n):
        a = xy[i]
        b = xy[(i + 1) % n]
        length = float(np.linalg.norm(b - a))
        nseg = max(1, int(np.round(length / h)))
        for k in range(nseg):
            parts.append(a + (b - a) * (k / nseg))
    return np.asarray(parts, float)


def _smooth_interior(cloud: np.ndarray, tris: np.ndarray, n_boundary: int,
                     xy: np.ndarray, h: float) -> np.ndarray | None:
    """Move interior vertices toward their neighbours, and keep them off the boundary.

    Each new position is the mean of the unique neighbours on the cloud before
    this pass, so the moves do not depend on each other.
    """
    if n_boundary >= len(cloud) or len(tris) == 0:
        return None
    tri = np.asarray(tris, np.int64)
    directed = np.concatenate([
        np.stack([tri[:, i], tri[:, j]], axis=1)
        for i, j in ((0, 1), (1, 2), (2, 0))
    ], axis=0)
    pairs = np.concatenate([directed, directed[:, ::-1]], axis=0)
    n = len(cloud)
    degree = np.bincount(pairs[:, 0], minlength=n)
    order = np.lexsort((pairs[:, 1], pairs[:, 0]))
    pairs = pairs[order]
    distinct = np.empty(len(pairs), dtype=bool)
    distinct[0] = True
    distinct[1:] = np.any(pairs[1:] != pairs[:-1], axis=1)
    pairs = pairs[distinct]
    src, dst = pairs[:, 0], pairs[:, 1]
    starts = np.concatenate([[0], np.flatnonzero(np.diff(src)) + 1])
    means = np.add.reduceat(cloud[dst], starts) / np.diff(np.append(starts, len(src)))[:, None]
    owners = src[starts]
    eligible = (owners >= n_boundary) & (degree[owners] >= 2)
    if not np.any(eligible):
        return None
    acc = means[eligible]
    margin = 0.35 * float(h)
    accept = points_in_ring(xy, acc) & (_edge_clearance(xy, acc) >= margin)
    if not np.any(accept):
        return None
    moved = np.asarray(cloud[n_boundary:], float).copy()
    moved[owners[eligible][accept] - n_boundary] = acc[accept]
    return moved


def _unstructured_seeds(xy: np.ndarray, h: float) -> np.ndarray:
    """Interior points about ``h`` apart, inset from the ring."""
    r = float(h)
    if r < 1e-15 or len(xy) < 3:
        return np.zeros((0, 2), float)
    margin = 0.62 * r
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    span = np.maximum(hi - lo, 1e-15)
    step = max(margin, 0.35 * r)
    xs = np.arange(float(lo[0]) + 0.5 * step, float(hi[0]), step)
    ys = np.arange(float(lo[1]) + 0.5 * step, float(hi[1]), step)
    if len(xs) == 0 or len(ys) == 0:
        return np.zeros((0, 2), float)
    xx, yy = np.meshgrid(xs, ys)
    cand = np.column_stack([xx.ravel(), yy.ravel()])
    rng = np.random.default_rng(20261117)
    cand = cand + rng.uniform(-0.48 * step, 0.48 * step, cand.shape)
    inside = points_in_ring(xy, cand)
    cand = cand[inside]
    if len(cand) == 0:
        return np.zeros((0, 2), float)
    clear = _edge_clearance(xy, cand) >= margin
    cand = cand[clear]
    if len(cand) == 0:
        return np.zeros((0, 2), float)
    order = rng.permutation(len(cand))
    # Grow from a candidate near the middle of the ring so the disk fills evenly.
    centre = xy.mean(axis=0)
    start_i = int(np.argmin(np.linalg.norm(cand - centre, axis=1)))
    order = np.concatenate([[start_i], order[order != start_i]])
    cell = r / np.sqrt(2.0)
    nx = max(int(np.ceil(span[0] / cell)) + 2, 1)
    ny = max(int(np.ceil(span[1] / cell)) + 2, 1)
    grid = -np.ones((ny, nx), np.int32)
    kept: list[np.ndarray] = []

    def slot(p: np.ndarray) -> tuple[int, int]:
        ix = int((p[0] - lo[0]) / cell)
        iy = int((p[1] - lo[1]) / cell)
        return int(np.clip(ix, 0, nx - 1)), int(np.clip(iy, 0, ny - 1))

    def room(p: np.ndarray) -> bool:
        ix, iy = slot(p)
        for yy in range(max(0, iy - 2), min(ny, iy + 3)):
            for xx in range(max(0, ix - 2), min(nx, ix + 3)):
                j = int(grid[yy, xx])
                if j >= 0 and float(np.linalg.norm(p - kept[j])) < r:
                    return False
        return True

    for i in order:
        p = cand[int(i)]
        if not room(p):
            continue
        ix, iy = slot(p)
        grid[iy, ix] = len(kept)
        kept.append(p)
    if not kept:
        return np.zeros((0, 2), float)
    return np.asarray(kept, float)


def _orient2(pts: np.ndarray, a: int, b: int, c: int) -> float:
    return float(
        (pts[b, 0] - pts[a, 0]) * (pts[c, 1] - pts[a, 1])
        - (pts[b, 1] - pts[a, 1]) * (pts[c, 0] - pts[a, 0])
    )


def _segments_cross(pts: np.ndarray, a: int, b: int, c: int, d: int) -> bool:
    """Proper intersection of open segments ab and cd."""
    o1 = _orient2(pts, a, b, c)
    o2 = _orient2(pts, a, b, d)
    o3 = _orient2(pts, c, d, a)
    o4 = _orient2(pts, c, d, b)
    return o1 * o2 < 0.0 and o3 * o4 < 0.0


def _quad_convex(pts: np.ndarray, a: int, b: int, c: int, d: int) -> bool:
    cycle = (a, b, c, d)
    signs = []
    for i in range(4):
        signs.append(_orient2(pts, cycle[i], cycle[(i + 1) % 4], cycle[(i + 2) % 4]))
    return all(s > 0.0 for s in signs) or all(s < 0.0 for s in signs)


def _recover_boundary(pts: np.ndarray, tris: np.ndarray, n_boundary: int) -> np.ndarray:
    """Flip until each consecutive boundary edge is a triangle edge.

    Boundary vertices stay as given. Interior points may sit in a boundary
    edge's diametral circle; the flip puts that edge back.
    """
    faces = [tuple(int(x) for x in t) for t in tris]
    if n_boundary < 3 or not faces:
        return tris

    def adjacency(faces: list[tuple[int, int, int]]) -> dict[tuple[int, int], list[int]]:
        edge_tri: dict[tuple[int, int], list[int]] = {}
        for ti, (a, b, c) in enumerate(faces):
            for u, v in ((a, b), (b, c), (c, a)):
                edge_tri.setdefault((u, v) if u < v else (v, u), []).append(ti)
        return edge_tri

    def positive(a: int, b: int, c: int) -> tuple[int, int, int]:
        if _orient2(pts, a, b, c) < 0.0:
            b, c = c, b
        return a, b, c

    edge_tri = adjacency(faces)
    for i in range(n_boundary):
        a = i
        b = (i + 1) % n_boundary
        key = (a, b) if a < b else (b, a)
        guard = 0
        limit = 8 * max(len(faces), 1)
        while key not in edge_tri and guard < limit:
            guard += 1
            crossed = None
            for (u, v), tids in edge_tri.items():
                if len(tids) != 2 or u in (a, b) or v in (a, b):
                    continue
                if not _segments_cross(pts, a, b, u, v):
                    continue
                t0, t1 = tids
                p = next(x for x in faces[t0] if x != u and x != v)
                q = next(x for x in faces[t1] if x != u and x != v)
                if _quad_convex(pts, u, p, v, q):
                    crossed = (u, v, p, q, t0, t1)
                    break
            if crossed is None:
                break
            u, v, p, q, t0, t1 = crossed
            faces[t0] = positive(p, u, q)
            faces[t1] = positive(p, q, v)
            edge_tri = adjacency(faces)
    return np.asarray(faces, np.int32)


def _triangle_array(parts: list[np.ndarray]) -> np.ndarray:
    if not parts:
        return np.zeros((0, 3, 3), float)
    return np.stack([np.asarray(part, float) for part in parts], axis=0)


def triangulate_ring(
    ring: np.ndarray, *, size: float | None = None, engine: str | None = None,
    pattern: str = "hex", split_boundary: bool = True,
) -> np.ndarray:
    """Delaunay triangles of a planar ring, one ``(n, 3, 3)`` array.

    ``pattern="hex"`` is the regular Steiner lattice. ``pattern="unstructured"``
    places interior points at ``size`` (or the ring's own edge length) and
    smooths them. With ``split_boundary`` the boundary edges are split to that
    spacing too; the original ring vertices stay either way. ``size`` is a
    minimum on the hex lattice; on the unstructured fill it is the spacing.
    """
    pts = np.asarray(ring, float)
    if len(pts) < 3:
        return np.zeros((0, 3, 3), float)
    origin, u, v, xy = _ring_basis(pts)
    xy = _unique_xy(xy)
    if len(xy) < 3:
        return np.zeros((0, 3, 3), float)
    if _area2(xy) < 0.0:
        pts = pts[::-1].copy()
        xy = xy[::-1].copy()
    if pattern == "hex":
        h = _target_size(xy, size)
        boundary = xy
        seeds = _hex_seeds(xy, h) if h >= 1e-12 else np.zeros((0, 2), float)
    elif pattern == "unstructured":
        h = float(size) if size is not None else _median_edge(xy)
        boundary = _split_boundary(xy, h) if split_boundary and h >= 1e-12 else xy
        seeds = _unstructured_seeds(boundary, h) if h >= 1e-12 else np.zeros((0, 2), float)
    else:
        raise ValueError(f"unknown triangulation pattern {pattern!r}")
    if h < 1e-12:
        return _triangle_array(_fan(pts))
    cloud = boundary if len(seeds) == 0 else np.vstack([boundary, seeds])
    try:
        tris = _simplices(cloud, engine)
    except Exception:
        return _triangle_array(_fan(pts))
    if len(tris) == 0:
        return _triangle_array(_fan(pts))
    if pattern == "unstructured":
        tris = _recover_boundary(cloud, tris, len(boundary))
        for _ in range(4):
            moved = _smooth_interior(cloud, tris, len(boundary), boundary, h)
            if moved is None:
                break
            seeds = moved
            cloud = np.vstack([boundary, seeds])
            try:
                tris = _simplices(cloud, engine)
            except Exception:
                break
            if len(tris) == 0:
                break
            tris = _recover_boundary(cloud, tris, len(boundary))
    cents = (cloud[tris[:, 0]] + cloud[tris[:, 1]] + cloud[tris[:, 2]]) / 3.0
    keep = points_in_ring(boundary, cents)
    a = cloud[tris[:, 0]]
    b = cloud[tris[:, 1]]
    c = cloud[tris[:, 2]]
    area = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    keep = keep & (np.abs(area) > 1e-16)
    idx = np.nonzero(keep)[0]
    if len(idx) == 0:
        return _triangle_array(_fan(pts))
    kept = np.stack([a[idx], b[idx], c[idx]], axis=1)
    return origin + kept[:, :, 0, None] * u + kept[:, :, 1, None] * v


def triangulate_rings(rings, *, size: float | None = None, engine: str | None = None,
                      pattern: str = "hex") -> list[np.ndarray]:
    return [triangulate_ring(ring, size=size, engine=engine, pattern=pattern) for ring in rings]


def fill_rib_cells(
    cells: list[tuple[str, np.ndarray]], *, engine: str | None = None,
) -> list[tuple[str, np.ndarray]]:
    """Triangle fill of already-extracted cell rings."""
    out: list[tuple[str, np.ndarray]] = []
    for name, ring in cells:
        for face in triangulate_ring(ring, engine=engine):
            out.append((name, face))
    return out


def mesh_rib_cells_many(
    template, samples: Mapping[str, np.ndarray], members, *, engine: str | None = None,
) -> list[tuple[str, np.ndarray]]:
    """Native triangle fill of every authored rib cell on ``members``."""
    from implicit_geometry_kernel.blade.topology.ribs import rib_cell_polygons_from_samples

    out: list[tuple[str, np.ndarray]] = []
    for member in members:
        for name, ring in rib_cell_polygons_from_samples(template, samples, member):
            for face in triangulate_ring(ring, engine=engine):
                out.append((name, face))
    return out


def mesh_rib_cells(template, samples: Mapping[str, np.ndarray], member, *, engine: str | None = None) -> list[tuple[str, np.ndarray]]:
    """Native triangle fill of each authored rib cell."""
    return mesh_rib_cells_many(template, samples, (member,), engine=engine)
