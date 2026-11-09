"""Unstructured polygonal fill of a planar n-gon (Voronoi clip, >4 nodes).

Kernel-owned mesher: no Gmsh. Element size follows the boundary sampling.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np
from numba import njit

from implicit_geometry_kernel.geos import lib as geos_lib
from implicit_geometry_kernel.geos.geom import MultiPolygon, Polygon
from implicit_geometry_kernel.geos.wkb import dumps, loads


def _ring_basis(ring: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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


def _area2(xy: np.ndarray) -> float:
    x, y = xy[:, 0], xy[:, 1]
    return float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _cross2(a, b, c) -> float:
    return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def points_in_ring(xy: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Winding test for many points against a planar ring (on-edge counts as in)."""
    xy = np.asarray(xy, float)
    pts = np.asarray(pts, float).reshape(-1, 2)
    if len(pts) == 0 or len(xy) < 3:
        return np.zeros(len(pts), dtype=bool)
    x1, y1 = xy[:, 0], xy[:, 1]
    x2, y2 = np.roll(x1, -1), np.roll(y1, -1)
    px = pts[:, 0][:, None]
    py = pts[:, 1][:, None]
    cross = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1)
    on = (
        (np.abs(cross) <= 1e-12)
        & (np.minimum(x1, x2) - 1e-12 <= px)
        & (px <= np.maximum(x1, x2) + 1e-12)
        & (np.minimum(y1, y2) - 1e-12 <= py)
        & (py <= np.maximum(y1, y2) + 1e-12)
    )
    hits = ((y1 > py) != (y2 > py)) & (px < (x2 - x1) * (py - y1) / (y2 - y1 + 1e-30) + x1)
    return on.any(axis=1) | (hits.sum(axis=1) % 2 == 1)


def point_in_ring(p, xy: np.ndarray) -> bool:
    return bool(points_in_ring(xy, np.asarray(p, float).reshape(1, 2))[0])


def _closed_xy(xy: np.ndarray) -> np.ndarray:
    xy = np.asarray(xy, float)
    if len(xy) and not np.allclose(xy[0], xy[-1], atol=1e-12, rtol=0.0):
        xy = np.vstack([xy, xy[0]])
    return xy


def _median_edge(xy: np.ndarray) -> float:
    d = np.linalg.norm(np.roll(xy, -1, axis=0) - xy, axis=1)
    d = d[d > 1e-15]
    return float(np.median(d)) if len(d) else 0.0


def _target_size(xy: np.ndarray, size: float | None) -> float:
    area = 0.5 * abs(_area2(xy))
    h = float(size) if size is not None else _median_edge(xy)
    if area > 0.0:
        h = max(h, float(np.sqrt(area / 96.0)))
    return h


def _hex_seeds(xy: np.ndarray, h: float) -> np.ndarray:
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    hx, hy = float(h), float(h) * 0.8660254037844386
    x0 = float(lo[0]) + 0.45 * h
    y0 = float(lo[1]) + 0.45 * h
    x1 = float(hi[0]) - 0.35 * h
    y1 = float(hi[1]) - 0.35 * h
    if hx < 1e-15 or hy < 1e-15 or x1 < x0 or y1 < y0:
        return np.zeros((0, 2), float)
    nx = int(np.floor((x1 - x0) / hx)) + 1
    ny = int(np.floor((y1 - y0) / hy)) + 1
    xs = x0 + np.arange(nx) * hx
    ys = y0 + np.arange(ny) * hy
    xx, yy = np.meshgrid(xs, ys)
    xx[1::2] += 0.5 * hx
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    pts = pts[(pts[:, 0] <= x1 + 1e-12) & (pts[:, 1] <= y1 + 1e-12)]
    if len(pts) == 0:
        return np.zeros((0, 2), float)
    return pts[points_in_ring(xy, pts)]


def _geom_rings(geom) -> list[np.ndarray]:
    if geom is None or getattr(geom, "is_empty", True):
        return []
    parts = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
    out: list[np.ndarray] = []
    for part in parts:
        if not isinstance(part, Polygon) or part.is_empty:
            continue
        xy = np.asarray(part.exterior.coords, float)[:, :2]
        if len(xy) > 2 and np.allclose(xy[0], xy[-1], atol=1e-12, rtol=0.0):
            xy = xy[:-1]
        if len(xy) >= 3:
            out.append(xy)
    return out


def _geos_clip(cell: np.ndarray, domain: Polygon) -> list[np.ndarray]:
    if len(cell) < 3:
        return []
    got = loads(geos_lib.intersection(dumps(Polygon(_closed_xy(cell))), dumps(domain)))
    return _geom_rings(got)


@njit(cache=True)
def _seg_hit(ax, ay, bx, by, cx, cy, dx, dy):
    den = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx)
    if abs(den) < 1e-18:
        return False, 0.0, 0.0
    t = ((cx - ax) * (dy - cy) - (cy - ay) * (dx - cx)) / den
    u = ((cx - ax) * (by - ay) - (cy - ay) * (bx - ax)) / den
    if t < -1e-12 or t > 1.0 + 1e-12 or u < -1e-12 or u > 1.0 + 1e-12:
        return False, 0.0, 0.0
    return True, ax + t * (bx - ax), ay + t * (by - ay)


@njit(cache=True)
def _convex_has(poly, n, px, py) -> bool:
    for i in range(n):
        ax, ay = poly[i, 0], poly[i, 1]
        bx, by = poly[(i + 1) % n, 0], poly[(i + 1) % n, 1]
        if (bx - ax) * (py - ay) - (by - ay) * (px - ax) < -1e-12:
            return False
    return True


@njit(cache=True)
def _ring_has(ring, n, px, py) -> bool:
    inside = False
    for i in range(n):
        x1, y1 = ring[i, 0], ring[i, 1]
        x2, y2 = ring[(i + 1) % n, 0], ring[(i + 1) % n, 1]
        cross = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1)
        if abs(cross) <= 1e-12:
            if min(x1, x2) - 1e-12 <= px <= max(x1, x2) + 1e-12 and (
                min(y1, y2) - 1e-12 <= py <= max(y1, y2) + 1e-12
            ):
                return True
        if ((y1 > py) != (y2 > py)) and (px < (x2 - x1) * (py - y1) / (y2 - y1 + 1e-30) + x1):
            inside = not inside
    return inside


@njit(cache=True)
def _clip_convex_to_ring(cell, ring, sx, sy):
    """Intersection of a convex Voronoi cell with a (possibly concave) ring."""
    nc = cell.shape[0]
    nr = ring.shape[0]
    buf = np.empty((nc + nr + nc * nr, 2), np.float64)
    k = 0
    for i in range(nc):
        if _ring_has(ring, nr, cell[i, 0], cell[i, 1]):
            buf[k, 0], buf[k, 1] = cell[i, 0], cell[i, 1]
            k += 1
    for j in range(nr):
        if _convex_has(cell, nc, ring[j, 0], ring[j, 1]):
            buf[k, 0], buf[k, 1] = ring[j, 0], ring[j, 1]
            k += 1
    for i in range(nc):
        ax, ay = cell[i, 0], cell[i, 1]
        bx, by = cell[(i + 1) % nc, 0], cell[(i + 1) % nc, 1]
        for j in range(nr):
            hit, x, y = _seg_hit(
                ax, ay, bx, by, ring[j, 0], ring[j, 1],
                ring[(j + 1) % nr, 0], ring[(j + 1) % nr, 1],
            )
            if hit:
                buf[k, 0], buf[k, 1] = x, y
                k += 1
    if k < 3:
        return np.empty((0, 2), np.float64)
    keep = np.ones(k, np.uint8)
    for i in range(k):
        if keep[i] == 0:
            continue
        for j in range(i + 1, k):
            dx = buf[i, 0] - buf[j, 0]
            dy = buf[i, 1] - buf[j, 1]
            if dx * dx + dy * dy < 1e-24:
                keep[j] = 0
    m = 0
    for i in range(k):
        if keep[i]:
            buf[m, 0], buf[m, 1] = buf[i, 0], buf[i, 1]
            m += 1
    k = m
    if k < 3:
        return np.empty((0, 2), np.float64)
    ang = np.empty(k, np.float64)
    for i in range(k):
        ang[i] = np.arctan2(buf[i, 1] - sy, buf[i, 0] - sx)
    for i in range(1, k):
        key_a = ang[i]
        key_x, key_y = buf[i, 0], buf[i, 1]
        j = i
        while j > 0 and ang[j - 1] > key_a:
            ang[j] = ang[j - 1]
            buf[j, 0], buf[j, 1] = buf[j - 1, 0], buf[j - 1, 1]
            j -= 1
        ang[j] = key_a
        buf[j, 0], buf[j, 1] = key_x, key_y
    return buf[:k].copy()


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
