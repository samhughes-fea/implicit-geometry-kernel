"""Gmsh Frontal-Delaunay fill of a planar n-gon."""

from __future__ import annotations

from typing import Mapping

import numpy as np


def _gmsh_start():
    import gmsh

    own = False
    if not gmsh.isInitialized():
        gmsh.initialize(["gmsh", "-v", "0"])
        own = True
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    return gmsh, own


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


def _fill_model(xy: np.ndarray, h: float) -> list[np.ndarray]:
    import gmsh

    n = len(xy)
    tags = [gmsh.model.geo.addPoint(float(p[0]), float(p[1]), 0.0, h) for p in xy]
    lines = [gmsh.model.geo.addLine(tags[i], tags[(i + 1) % n]) for i in range(n)]
    loop = gmsh.model.geo.addCurveLoop(lines)
    gmsh.model.geo.addPlaneSurface([loop])
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    node_tags, coord, _ = gmsh.model.mesh.getNodes()
    xyz = np.asarray(coord, float).reshape(-1, 3)
    lookup = {int(t): xyz[i] for i, t in enumerate(node_tags)}
    types, _etags, ntags = gmsh.model.mesh.getElements(2)
    tris: list[np.ndarray] = []
    for typ, nodes in zip(types, ntags):
        if int(typ) != 2:
            continue
        arr = np.asarray(nodes, int).reshape(-1, 3)
        for row in arr:
            tris.append(np.vstack([lookup[int(t)][:2] for t in row]))
    return tris


def triangulate_ring(ring: np.ndarray, *, size: float | None = None) -> list[np.ndarray]:
    """Frontal-Delaunay triangles of a planar ring; size defaults to median edge."""
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
    gmsh, own = _gmsh_start()
    try:
        gmsh.model.add("rib_cell")
        xy_tris = _fill_model(xy, h)
    except Exception:
        return _fan(pts)
    finally:
        if own:
            gmsh.finalize()
    out = []
    for cell in xy_tris:
        if abs(_area2(cell)) < 1e-16:
            continue
        out.append(origin + np.outer(cell[:, 0], u) + np.outer(cell[:, 1], v))
    return out or _fan(pts)


def triangulate_rings(rings, *, size: float | None = None) -> list[list[np.ndarray]]:
    """One Gmsh session for many planar rings."""
    gmsh, own = _gmsh_start()
    try:
        out: list[list[np.ndarray]] = []
        for i, ring in enumerate(rings):
            gmsh.model.add(f"rib_{i}")
            gmsh.model.setCurrent(f"rib_{i}")
            pts = np.asarray(ring, float)
            if len(pts) < 3:
                out.append([])
                continue
            origin, u, v, xy = _ring_basis(pts)
            xy = _unique_xy(xy)
            if len(xy) < 3:
                out.append([])
                continue
            if _area2(xy) < 0.0:
                pts = pts[::-1].copy()
                xy = xy[::-1].copy()
            h = _target_size(xy, size)
            if h < 1e-12:
                out.append([pts.copy()])
                continue
            try:
                xy_tris = _fill_model(xy, h)
            except Exception:
                xy_tris = None
            try:
                gmsh.model.remove()
            except Exception:
                pass
            if xy_tris is None:
                out.append(_fan(pts))
                continue
            faces = []
            for cell in xy_tris:
                if abs(_area2(cell)) < 1e-16:
                    continue
                faces.append(origin + np.outer(cell[:, 0], u) + np.outer(cell[:, 1], v))
            out.append(faces or _fan(pts))
        return out
    finally:
        if own:
            gmsh.finalize()


def fill_rib_cells(cells: list[tuple[str, np.ndarray]]) -> list[tuple[str, np.ndarray]]:
    """One Gmsh session for already-extracted cell rings."""
    if not cells:
        return []
    faces = triangulate_rings([ring for _name, ring in cells])
    out: list[tuple[str, np.ndarray]] = []
    for (name, _ring), tris in zip(cells, faces):
        for face in tris:
            out.append((name, face))
    return out


def mesh_rib_cells_many(template, samples: Mapping[str, np.ndarray], members) -> list[tuple[str, np.ndarray]]:
    """One Gmsh session for every authored rib cell on ``members``."""
    from implicit_geometry_kernel.blade.topology.ribs import rib_cell_polygons_from_samples

    cells: list[tuple[str, np.ndarray]] = []
    for member in members:
        cells.extend(rib_cell_polygons_from_samples(template, samples, member))
    return fill_rib_cells(cells)


def mesh_rib_cells(template, samples: Mapping[str, np.ndarray], member) -> list[tuple[str, np.ndarray]]:
    """Gmsh triangle fill of each authored rib cell (own session)."""
    return mesh_rib_cells_many(template, samples, (member,))
