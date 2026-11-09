"""Steiner hex + Delaunay fill of a planar n-gon (Qhull, else numba)."""

from __future__ import annotations

from typing import Mapping

import numpy as np
from numba import njit

from implicit_geometry_kernel.mesh.polygon import _hex_seeds, _ring_basis, points_in_ring


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


def triangulate_ring(
    ring: np.ndarray, *, size: float | None = None, engine: str | None = None,
) -> list[np.ndarray]:
    """Delaunay triangles of a planar ring with hex Steiner points."""
    pts = np.asarray(ring, float)
    if len(pts) < 3:
        return []
    origin, u, v, xy = _ring_basis(pts)
    xy = _unique_xy(xy)
    if len(xy) < 3:
        return []
    if _area2(xy) < 0.0:
        pts = pts[::-1].copy()
        xy = xy[::-1].copy()
    h = _target_size(xy, size)
    if h < 1e-12:
        return [pts.copy()]
    seeds = _hex_seeds(xy, h)
    cloud = xy if len(seeds) == 0 else np.vstack([xy, seeds])
    try:
        tris = _simplices(cloud, engine)
    except Exception:
        return _fan(pts)
    if len(tris) == 0:
        return _fan(pts)
    cents = (cloud[tris[:, 0]] + cloud[tris[:, 1]] + cloud[tris[:, 2]]) / 3.0
    keep = points_in_ring(xy, cents)
    a = cloud[tris[:, 0]]
    b = cloud[tris[:, 1]]
    c = cloud[tris[:, 2]]
    area = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    keep = keep & (np.abs(area) > 1e-16)
    idx = np.nonzero(keep)[0]
    if len(idx) == 0:
        return _fan(pts)
    kept = np.stack([a[idx], b[idx], c[idx]], axis=1)
    xyz = origin + kept[:, :, 0, None] * u + kept[:, :, 1, None] * v
    return [xyz[i] for i in range(len(xyz))]


def triangulate_rings(rings, *, size: float | None = None, engine: str | None = None) -> list[list[np.ndarray]]:
    return [triangulate_ring(ring, size=size, engine=engine) for ring in rings]


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
