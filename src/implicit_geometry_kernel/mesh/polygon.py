"""Unstructured polygonal fill of a planar n-gon (conforming Voronoi, >4 nodes).

Kernel-owned mesher: no Gmsh. Element size follows the boundary sampling.
"""

from __future__ import annotations


import numpy as np

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



def polygon_mesh_ring(*args, **kwargs):
    from implicit_geometry_kernel.mesh.polygon_voronoi import polygon_mesh_ring as _fn
    return _fn(*args, **kwargs)

def fill_rib_cells(*args, **kwargs):
    from implicit_geometry_kernel.mesh.polygon_voronoi import fill_rib_cells as _fn
    return _fn(*args, **kwargs)

def mesh_rib_cells(*args, **kwargs):
    from implicit_geometry_kernel.mesh.polygon_voronoi import mesh_rib_cells as _fn
    return _fn(*args, **kwargs)

def mesh_rib_cells_many(*args, **kwargs):
    from implicit_geometry_kernel.mesh.polygon_voronoi import mesh_rib_cells_many as _fn
    return _fn(*args, **kwargs)
