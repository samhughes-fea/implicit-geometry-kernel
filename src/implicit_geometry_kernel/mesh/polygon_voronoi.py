"""Voronoi cells and n-gon fill (Qhull or Bowyer–Watson dual)."""

from __future__ import annotations

from typing import Mapping

import numpy as np

from implicit_geometry_kernel.geos.geom import Polygon
from implicit_geometry_kernel.mesh.polygon import (
    _area2,
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


def _clip_cells(raw: list[np.ndarray], ring: np.ndarray, domain: Polygon) -> list[np.ndarray]:
    """Exact GEOS intersection of each Voronoi cell with the ring (CCW pieces).

    A cell cut by a concave ring can fall in more than one piece; each piece is
    its own cell. An angle sort about the seed is only right for star-shaped
    pieces, which thin trailing-edge cells are not.
    """
    out: list[np.ndarray] = []
    for cell in raw:
        if len(cell) < 3:
            continue
        if bool(np.all(points_in_ring(ring, cell))):
            pieces = [cell]
            if not _inside_exactly(cell, domain):
                pieces = _geos_clip(cell, domain)
        else:
            pieces = _geos_clip(cell, domain)
        for piece in pieces:
            piece = np.asarray(piece, float)
            if len(piece) >= 3 and abs(_area2(piece)) > 1e-300:
                out.append(piece if _area2(piece) > 0.0 else piece[::-1].copy())
    return out


def _inside_exactly(cell: np.ndarray, domain: Polygon) -> bool:
    """All vertices inside is not enough: a concave ring can cut an edge."""
    got = _geos_clip(cell, domain)
    return len(got) == 1 and abs(abs(_area2(got[0])) - abs(_area2(cell))) <= 1e-12 * abs(_area2(cell))


def _ring_hit(p: np.ndarray, a: np.ndarray, ab: np.ndarray, l2: np.ndarray, tol: float) -> tuple[int, float] | None:
    t = np.clip(((p - a) * ab).sum(axis=1) / l2, 0.0, 1.0)
    d = np.linalg.norm(a + t[:, None] * ab - p, axis=1)
    i = int(np.argmin(d))
    return (i, float(t[i])) if d[i] <= tol else None


def _snap_to_ring(cells: list[np.ndarray], ring: np.ndarray) -> tuple[np.ndarray, list[list[int]]]:
    """Cells as vertex indices whose boundary nodes are exactly the ring nodes.

    A clipped cell meets the ring at points along ring edges. Each such point
    moves to the nearer ring node, and the cell takes every ring node between
    its entry and exit, so a cell and the skin beside it share nodes and edges:
    no hanging node on the ring. Neighbouring cells snap a shared point the
    same way, so they stay conforming. Indices ``0 .. n-1`` are the ring nodes.
    """
    n = len(ring)
    a = ring
    ab = np.roll(ring, -1, axis=0) - ring
    l2 = np.maximum((ab * ab).sum(axis=1), 1e-300)
    tol = 1e-9 * max(float(np.sqrt(l2.max())), 1e-300)
    verts: list[np.ndarray] = [p for p in ring]
    index: dict[tuple[float, float], int] = {}

    def interior(p: np.ndarray) -> int:
        key = (round(float(p[0]), 12), round(float(p[1]), 12))
        if key not in index:
            index[key] = len(verts)
            verts.append(np.asarray(p, float))
        return index[key]

    faces: list[list[int]] = []
    for cell in cells:
        m = len(cell)
        tokens: list[int] = []
        on: list[bool] = []
        for p in cell:
            hit = _ring_hit(p, a, ab, l2, tol)
            if hit is None:
                tokens.append(interior(p))
                on.append(False)
            else:
                tokens.append(hit[0] if hit[1] < 0.5 else (hit[0] + 1) % n)
                on.append(True)
        seq: list[int] = []
        for k in range(m):
            seq.append(tokens[k])
            j = (k + 1) % m
            if on[k] and on[j] and tokens[k] != tokens[j]:
                mid = 0.5 * (cell[k] + cell[j])
                if _ring_hit(mid, a, ab, l2, tol) is not None:
                    r = (tokens[k] + 1) % n
                    while r != tokens[j]:
                        seq.append(r)
                        r = (r + 1) % n
        face: list[int] = []
        for v in seq:
            if not face or face[-1] != v:
                face.append(v)
        while len(face) > 1 and face[0] == face[-1]:
            face.pop()
        if len(face) >= 3:
            faces.append(face)
    return np.asarray(verts, float), faces


def _face_area(verts: np.ndarray, face: list[int]) -> float:
    return 0.5 * _area2(verts[face])


def _tiles(verts: np.ndarray, faces: list[list[int]], ring_area: float) -> bool:
    """Positive cells whose areas add up to the ring: an exact, unfolded tiling."""
    if not faces:
        return False
    areas = np.array([_face_area(verts, f) for f in faces])
    return bool(np.all(areas > 0.0)) and abs(float(areas.sum()) - ring_area) <= 1e-9 * ring_area


def _merge_edge_neighbour(faces: list[list[int]], k: int, verts: np.ndarray) -> bool:
    """Merge face ``k`` into the neighbour it shares the longest edge run with."""
    edges_k = {(face_a, face_b) for face_a, face_b in zip(faces[k], faces[k][1:] + faces[k][:1])}
    best, best_len = -1, 0.0
    for j, other in enumerate(faces):
        if j == k:
            continue
        shared = [(b, a) for a, b in zip(other, other[1:] + other[:1]) if (b, a) in edges_k]
        if not shared:
            continue
        length = sum(float(np.linalg.norm(verts[a] - verts[b])) for a, b in shared)
        if length > best_len:
            best, best_len = j, length
    if best < 0:
        return False
    merged = _union_walk(faces[k], faces[best])
    if merged is None or _face_area(verts, merged) <= 0.0:
        return False
    keep = [f for i, f in enumerate(faces) if i not in (k, best)]
    faces[:] = keep + [merged]
    return True


def _union_walk(f: list[int], g: list[int]) -> list[int] | None:
    """Boundary of two CCW faces that share a run of edges."""
    edges = list(zip(f, f[1:] + f[:1])) + list(zip(g, g[1:] + g[:1]))
    es = set(edges)
    outer = [(a, b) for a, b in edges if (b, a) not in es]
    nxt: dict[int, int] = {}
    for a, b in outer:
        if a in nxt:
            return None
        nxt[a] = b
    if not outer:
        return None
    start = outer[0][0]
    walk = [start]
    v = nxt[start]
    while v != start:
        if v not in nxt or len(walk) > len(outer):
            return None
        walk.append(v)
        v = nxt[v]
    return walk if len(walk) == len(outer) else None


def _absorb_small(faces: list[list[int]], verts: np.ndarray) -> None:
    """3- and 4-node cells join an edge neighbour, so every cell has more than 4 nodes."""
    changed = True
    while changed:
        changed = False
        for k, face in enumerate(faces):
            if len(face) <= 4 and len(faces) > 1 and _merge_edge_neighbour(faces, k, verts):
                changed = True
                break


def _conforming_cells(xy: np.ndarray, h: float, lloyd: int, domain: Polygon) -> tuple[np.ndarray, list[list[int]]] | None:
    seeds = _hex_seeds(xy, h)
    if len(seeds) < 3:
        return None
    cells = _clip_cells(_voronoi_cells(seeds), xy, domain)
    for _ in range(max(int(lloyd), 0)):
        nxt = [c for c in (_cell_centroid(cell) for cell in cells) if point_in_ring(c, xy)]
        if len(nxt) < 3:
            break
        cells = _clip_cells(_voronoi_cells(np.asarray(nxt, float)), xy, domain)
    verts, faces = _snap_to_ring(cells, xy)
    if not _tiles(verts, faces, 0.5 * _area2(xy)):
        return None
    _absorb_small(faces, verts)
    return verts, faces


def polygon_mesh_ring(ring: np.ndarray, *, size: float | None = None, lloyd: int = 0) -> list[np.ndarray]:
    """Fill a planar 3-D n-gon with unstructured n-gons (more than 4 nodes).

    The fill conforms to the ring: its boundary nodes are exactly the ring's
    nodes (returned bit for bit), and every ring edge is an edge of one cell.
    Interior cells are clipped Voronoi cells of a hexagonal seed grid of
    spacing ``size`` (default: the median ring edge). A conforming cell cannot
    be finer along the ring than the ring itself; if snapping to the ring would
    fold a cell, the spacing grows by 1.25 until the cells tile the ring.
    """
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
    whole = [pts.copy()] if len(pts) > 4 else []
    if h < 1e-12 or area < 1e-18:
        return whole
    span = float(np.linalg.norm(xy.max(axis=0) - xy.min(axis=0)))
    got = None
    while got is None and h <= span:
        got = _conforming_cells(xy, h, lloyd, domain)
        h *= 1.25
    if got is None:
        return whole
    verts, faces = got
    n = len(pts)
    out = []
    for face in faces:
        if len(face) <= 4 and len(faces) > 1:
            continue
        rows = [pts[i] if i < n else origin + verts[i, 0] * u + verts[i, 1] * v for i in face]
        out.append(np.asarray(rows, float))
    return out or whole


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
