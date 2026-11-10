"""Voronoi cells and n-gon fill (Qhull or Bowyer–Watson dual)."""

from __future__ import annotations

from typing import Mapping

import numpy as np

from implicit_geometry_kernel.geos.geom import Polygon
from implicit_geometry_kernel.mesh.polygon import (
    _area2,
    _clip_convex_to_ring,
    _closed_xy,
    _geos_clip,
    _hex_seeds,
    _ring_basis,
    _target_size,
    point_in_ring,
    points_in_ring,
)


def _has_scipy() -> bool:
    import importlib.util

    return importlib.util.find_spec("scipy") is not None


def _voronoi_cells_qhull(seeds: np.ndarray) -> list[np.ndarray]:
    from scipy.spatial import Voronoi

    lo, hi = seeds.min(axis=0), seeds.max(axis=0)
    d = 40.0 * (float(np.linalg.norm(hi - lo)) + 1.0)
    mid = 0.5 * (lo + hi)
    dummy = np.array(
        [[mid[0], mid[1] + d], [mid[0] - d, mid[1] - d], [mid[0] + d, mid[1] - d],
         [mid[0] + d, mid[1] + d], [mid[0] - d, mid[1] + d]],
        float,
    )
    vor = Voronoi(np.vstack([seeds, dummy]))
    out: list[np.ndarray] = []
    for i in range(len(seeds)):
        reg = vor.regions[int(vor.point_region[i])]
        if not reg or -1 in reg:
            out.append(np.zeros((0, 2), float))
            continue
        verts = np.asarray(vor.vertices[reg], float)
        if _area2(verts) < 0.0:
            verts = verts[::-1].copy()
        out.append(verts)
    return out


def _circumcenter(a, b, c) -> np.ndarray:
    d = 2.0 * ((a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])))
    if abs(d) < 1e-18:
        return (a + b + c) / 3.0
    a2, b2, c2 = np.dot(a, a), np.dot(b, b), np.dot(c, c)
    ux = (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d
    uy = (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d
    return np.array([ux, uy], float)


def _voronoi_cells_dual(seeds: np.ndarray) -> list[np.ndarray]:
    from implicit_geometry_kernel.mesh.delaunay import _bowyer_watson_numba

    tris = _bowyer_watson_numba(np.asarray(seeds, np.float64))
    n = len(seeds)
    buckets: list[list[np.ndarray]] = [[] for _ in range(n)]
    for t in tris:
        cc = _circumcenter(seeds[t[0]], seeds[t[1]], seeds[t[2]])
        for i in t:
            buckets[int(i)].append(cc)
    out: list[np.ndarray] = []
    for i, pts in enumerate(buckets):
        if len(pts) < 3:
            out.append(np.zeros((0, 2), float))
            continue
        arr = np.asarray(pts, float)
        ang = np.arctan2(arr[:, 1] - seeds[i, 1], arr[:, 0] - seeds[i, 0])
        arr = arr[np.argsort(ang)]
        if _area2(arr) < 0.0:
            arr = arr[::-1].copy()
        out.append(arr)
    return out


def _voronoi_cells(seeds: np.ndarray) -> list[np.ndarray]:
    if _has_scipy():
        try:
            return _voronoi_cells_qhull(seeds)
        except Exception:
            pass
    return _voronoi_cells_dual(seeds)


def _cell_centroid(xy: np.ndarray) -> np.ndarray:
    a2 = _area2(xy)
    if abs(a2) < 1e-18:
        return xy.mean(axis=0)
    x, y = xy[:, 0], xy[:, 1]
    cx = np.sum((x + np.roll(x, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)) / (3.0 * a2)
    cy = np.sum((y + np.roll(y, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)) / (3.0 * a2)
    return np.array([float(cx), float(cy)], float)


def _clip_one(cell: np.ndarray, ring: np.ndarray, seed: np.ndarray, domain: Polygon) -> np.ndarray:
    if len(cell) < 3:
        return np.zeros((0, 2), float)
    inside = points_in_ring(ring, cell)
    if bool(np.all(inside)):
        return cell
    got = _clip_convex_to_ring(np.ascontiguousarray(cell, np.float64),
                               np.ascontiguousarray(ring, np.float64),
                               float(seed[0]), float(seed[1]))
    if len(got) >= 3:
        return got
    repaired = _geos_clip(cell, domain)
    return repaired[0] if repaired else np.zeros((0, 2), float)


def _merge_small(cells: list[np.ndarray], seeds: np.ndarray, ring: np.ndarray) -> list[np.ndarray]:
    """Absorb 3-gons and 4-gons into the nearest face with more than 4 nodes."""
    big: list[tuple[np.ndarray, np.ndarray]] = []
    small: list[np.ndarray] = []
    for cell, seed in zip(cells, seeds):
        if len(cell) > 4 and abs(_area2(cell)) > 1e-16:
            big.append((cell, seed))
        elif len(cell) >= 3 and abs(_area2(cell)) > 1e-16:
            small.append(cell)
    if not big:
        return [c for c in cells if len(c) > 4]
    for sm in small:
        c = _cell_centroid(sm)
        j = int(np.argmin([np.linalg.norm(_cell_centroid(b[0]) - c) for b in big]))
        cell, seed = big[j]
        merged = np.vstack([cell, sm])
        ang = np.arctan2(merged[:, 1] - seed[1], merged[:, 0] - seed[0])
        order = np.argsort(ang)
        uniq = [merged[order[0]]]
        for i in order[1:]:
            if float(np.linalg.norm(merged[i] - uniq[-1])) > 1e-12:
                uniq.append(merged[i])
        arr = np.asarray(uniq, float)
        if len(arr) > 2 and float(np.linalg.norm(arr[0] - arr[-1])) < 1e-12:
            arr = arr[:-1]
        clipped = _clip_convex_to_ring(
            np.ascontiguousarray(arr, np.float64),
            np.ascontiguousarray(ring, np.float64),
            float(seed[0]), float(seed[1]),
        )
        if len(clipped) > 4:
            big[j] = (clipped, seed)
    return [c for c, _s in big if len(c) > 4]


def polygon_mesh_ring(ring: np.ndarray, *, size: float | None = None, lloyd: int = 0) -> list[np.ndarray]:
    """Fill a planar 3-D n-gon with unstructured n-gons (more than 4 nodes)."""
    pts = np.asarray(ring, float)
    if len(pts) > 2 and np.allclose(pts[0], pts[-1], atol=1e-12, rtol=0.0):
        pts = pts[:-1]
    if len(pts) < 3:
        return []
    origin, u, v, xy = _ring_basis(pts)
    if _area2(xy) < 0.0:
        pts = pts[::-1].copy()
        xy = xy[::-1].copy()
    domain = Polygon(_closed_xy(xy))
    area = 0.5 * abs(_area2(xy))
    h = _target_size(xy, size)
    if h < 1e-12 or area < 1e-18:
        return [pts.copy()] if len(pts) > 4 else []
    seeds = _hex_seeds(xy, h)
    if len(seeds) < 3:
        return [pts.copy()] if len(pts) > 4 else []
    raw = _voronoi_cells(seeds)
    cells = [_clip_one(cell, xy, seeds[i], domain) for i, cell in enumerate(raw)]
    for _ in range(max(int(lloyd), 0)):
        nxt = []
        for cell in cells:
            if len(cell) < 3:
                continue
            c = _cell_centroid(cell)
            if point_in_ring(c, xy):
                nxt.append(c)
        if len(nxt) < 3:
            break
        seeds = np.asarray(nxt, float)
        raw = _voronoi_cells(seeds)
        cells = [_clip_one(cell, xy, seeds[i], domain) for i, cell in enumerate(raw)]
    cells = _merge_small(cells, seeds, xy)
    out = []
    for cell in cells:
        if len(cell) <= 4 or abs(_area2(cell)) < 1e-16:
            continue
        out.append(origin + np.outer(cell[:, 0], u) + np.outer(cell[:, 1], v))
    return out or ([pts.copy()] if len(pts) > 4 else [])


def fill_rib_cells(cells: list[tuple[str, np.ndarray]]) -> list[tuple[str, np.ndarray]]:
    """N-gon fill of already-extracted cell rings."""
    out: list[tuple[str, np.ndarray]] = []
    for name, ring in cells:
        for face in polygon_mesh_ring(ring):
            out.append((name, face))
    return out


def mesh_rib_cells_many(template, samples: Mapping[str, np.ndarray], members) -> list[tuple[str, np.ndarray]]:
    from implicit_geometry_kernel.blade.topology.ribs import rib_cell_polygons_from_samples

    cells: list[tuple[str, np.ndarray]] = []
    for member in members:
        cells.extend(rib_cell_polygons_from_samples(template, samples, member))
    return fill_rib_cells(cells)


def mesh_rib_cells(template, samples: Mapping[str, np.ndarray], member) -> list[tuple[str, np.ndarray]]:
    return mesh_rib_cells_many(template, samples, (member,))
