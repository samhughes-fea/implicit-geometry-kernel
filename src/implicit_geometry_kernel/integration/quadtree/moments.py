"""Exact area moments of rectangles, polygons and cut leaf cells."""

from __future__ import annotations

from typing import Callable

import numpy as np
from implicit_geometry_kernel.integration.quadtree.marching import _interpolate_zero_crossing
from implicit_geometry_kernel.integration.quadtree.tree import QuadCell


def _rectangle_moments(x0: float, x1: float, y0: float, y1: float) -> tuple[float, ...]:
    """Exact integrals over an axis-aligned rectangle (area and raw moments about origin)."""
    dx = x1 - x0
    dy = y1 - y0
    area = dx * dy
    sx = 0.5 * area * (x0 + x1)
    sy = 0.5 * area * (y0 + y1)
    sxx = area * (x0 * x0 + x0 * x1 + x1 * x1) / 3.0
    syy = area * (y0 * y0 + y0 * y1 + y1 * y1) / 3.0
    sxy = area * (x0 * y0 + 0.5 * (x0 * y1 + x1 * y0) + x1 * y1) / 4.0
    return area, sx, sy, sxx, syy, sxy


def _rectangle_moments_batch(
    x0: np.ndarray,
    x1: np.ndarray,
    y0: np.ndarray,
    y1: np.ndarray,
) -> np.ndarray:
    """Batched rectangle moments; returns ``(n, 6)`` raw moment rows."""
    dx = x1 - x0
    dy = y1 - y0
    area = dx * dy
    sx = 0.5 * area * (x0 + x1)
    sy = 0.5 * area * (y0 + y1)
    sxx = area * (x0 * x0 + x0 * x1 + x1 * x1) / 3.0
    syy = area * (y0 * y0 + y0 * y1 + y1 * y1) / 3.0
    sxy = area * (x0 * y0 + 0.5 * (x0 * y1 + x1 * y0) + x1 * y1) / 4.0
    return np.column_stack([area, sx, sy, sxx, syy, sxy])


def _polygon_moments(verts: np.ndarray) -> tuple[float, ...]:
    """Shoelace area and raw second moments about the origin for a simple polygon."""
    v = np.asarray(verts, dtype=float).reshape(-1, 2)
    if len(v) < 3:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    x = v[:, 0]
    y = v[:, 1]
    x1 = np.roll(x, -1)
    y1 = np.roll(y, -1)
    cross = x * y1 - x1 * y
    area = 0.5 * float(np.sum(cross))
    if abs(area) < 1e-30:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    sx = float(np.sum((x + x1) * cross) / 6.0)
    sy = float(np.sum((y + y1) * cross) / 6.0)
    sxx = float(np.sum((x * x + x * x1 + x1 * x1) * cross) / 12.0)
    syy = float(np.sum((y * y + y * y1 + y1 * y1) * cross) / 12.0)
    sxy = float(
        np.sum((x * y1 + 2.0 * x * y + 2.0 * x1 * y1 + x1 * y) * cross) / 24.0
    )
    return area, sx, sy, sxx, syy, sxy


def _inside_polygon_from_cell(cell: QuadCell, phi_corners: np.ndarray) -> np.ndarray | None:
    """Build the inside (phi<0) polygon for a boundary leaf via marching squares."""
    corners = cell.corners()
    phi = np.asarray(phi_corners, dtype=float).reshape(4)
    inside_mask = phi <= 0.0
    if inside_mask.all():
        return corners
    if not inside_mask.any():
        return None
    pts: list[np.ndarray] = []
    for i in range(4):
        if inside_mask[i]:
            pts.append(corners[i])
        j = (i + 1) % 4
        if inside_mask[i] != inside_mask[j]:
            pts.append(_interpolate_zero_crossing(corners[i], corners[j], phi[i], phi[j]))
    if len(pts) < 3:
        return None
    return np.asarray(pts, dtype=float)


def _child_cells(cell: QuadCell) -> list[QuadCell]:
    xm = 0.5 * (cell.x_min + cell.x_max)
    ym = 0.5 * (cell.y_min + cell.y_max)
    return [
        QuadCell(cell.x_min, xm, cell.y_min, ym, level=cell.level + 1),
        QuadCell(xm, cell.x_max, cell.y_min, ym, level=cell.level + 1),
        QuadCell(xm, cell.x_max, ym, cell.y_max, level=cell.level + 1),
        QuadCell(cell.x_min, xm, ym, cell.y_max, level=cell.level + 1),
    ]


def _accumulate_cell(
    acc: list[float],
    cell: QuadCell,
    phi_samples: np.ndarray,
    *,
    sampler: Callable[[np.ndarray], np.ndarray] | None = None,
    subdiv_depth: int = 0,
    max_subdiv: int = 2,
) -> None:
    phi = np.asarray(phi_samples, dtype=float).reshape(-1)
    corners_phi = phi[:4]
    center_phi = float(phi[4]) if len(phi) >= 5 else float(np.mean(corners_phi))

    if np.all(corners_phi > 0.0) and center_phi > 0.0:
        return
    if np.all(corners_phi <= 0.0) and center_phi <= 0.0:
        m = _rectangle_moments(cell.x_min, cell.x_max, cell.y_min, cell.y_max)
    elif np.all(corners_phi > 0.0) and center_phi <= 0.0 and sampler is not None and subdiv_depth < max_subdiv:
        for child in _child_cells(cell):
            pts = child.corners()
            ctr = np.array([[0.5 * (child.x_min + child.x_max), 0.5 * (child.y_min + child.y_max)]])
            child_phi = np.concatenate([np.asarray(sampler(pts), dtype=float).reshape(4), sampler(ctr)])
            _accumulate_cell(
                acc,
                child,
                child_phi,
                sampler=sampler,
                subdiv_depth=subdiv_depth + 1,
                max_subdiv=max_subdiv,
            )
        return
    else:
        poly = _inside_polygon_from_cell(cell, corners_phi)
        if poly is None:
            return
        m = _polygon_moments(poly)
    for i in range(6):
        acc[i] += m[i]


def _boundary_recipe(code: int) -> list[tuple[int, int]]:
    """Ordered polygon vertex recipe for a 4-bit inside pattern.

    Mirrors ``_inside_polygon_from_cell``: walk corners 0..3, emit each inside
    corner (kind 0) then the zero-crossing on edge ``(i, i+1)`` (kind 1) when the
    edge straddles the boundary.
    """
    inside = [(code >> i) & 1 for i in range(4)]
    recipe: list[tuple[int, int]] = []
    for i in range(4):
        if inside[i]:
            recipe.append((0, i))
        if inside[i] != inside[(i + 1) % 4]:
            recipe.append((1, i))
    return recipe


_BOUNDARY_RECIPES = {code: _boundary_recipe(code) for code in range(16)}


def _shoelace_moments_batch(verts: np.ndarray) -> np.ndarray:
    """Summed raw moments ``(6,)`` over a batch of equal-length polygons.

    ``verts`` has shape ``(m, L, 2)``. Uses the same shoelace formulas as
    :func:`_polygon_moments`, including its near-zero-area guard.
    """
    x = verts[:, :, 0]
    y = verts[:, :, 1]
    xn = np.roll(x, -1, axis=1)
    yn = np.roll(y, -1, axis=1)
    cross = x * yn - xn * y
    area = 0.5 * cross.sum(axis=1)
    sx = ((x + xn) * cross).sum(axis=1) / 6.0
    sy = ((y + yn) * cross).sum(axis=1) / 6.0
    sxx = ((x * x + x * xn + xn * xn) * cross).sum(axis=1) / 12.0
    syy = ((y * y + y * yn + yn * yn) * cross).sum(axis=1) / 12.0
    sxy = ((x * yn + 2.0 * x * y + 2.0 * xn * yn + xn * y) * cross).sum(axis=1) / 24.0
    moments = np.column_stack([area, sx, sy, sxx, syy, sxy])
    moments[np.abs(area) < 1e-30] = 0.0
    return moments.sum(axis=0)


def _boundary_moments_batch(bounds: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Accumulated raw moments ``(6,)`` for boundary leaves, grouped by pattern.

    Vectorized equivalent of looping ``_inside_polygon_from_cell`` +
    ``_polygon_moments`` over each boundary cell. ``bounds`` is ``(n, 4)``
    (x_min, x_max, y_min, y_max, matching ``QuadCell.corners()`` order) and
    ``phi`` is the ``(n, 4)`` corner field.
    """
    acc = np.zeros(6, dtype=float)
    n = bounds.shape[0]
    if n == 0:
        return acc
    x0, x1, y0, y1 = bounds[:, 0], bounds[:, 1], bounds[:, 2], bounds[:, 3]
    corners = [
        np.stack([x0, y0], axis=1),
        np.stack([x1, y0], axis=1),
        np.stack([x1, y1], axis=1),
        np.stack([x0, y1], axis=1),
    ]
    inside = phi <= 0.0
    codes = (
        inside[:, 0].astype(np.int64)
        | (inside[:, 1].astype(np.int64) << 1)
        | (inside[:, 2].astype(np.int64) << 2)
        | (inside[:, 3].astype(np.int64) << 3)
    )
    edges: list[np.ndarray] = []
    for i in range(4):
        j = (i + 1) % 4
        va, vb = phi[:, i], phi[:, j]
        denom = vb - va
        degenerate = np.abs(denom) < 1e-15
        t = np.where(degenerate, 0.5, -va / np.where(degenerate, 1.0, denom))
        edges.append(corners[i] + t[:, None] * (corners[j] - corners[i]))
    for code in range(16):
        recipe = _BOUNDARY_RECIPES[code]
        if len(recipe) < 3:
            continue
        sel = np.flatnonzero(codes == code)
        if sel.size == 0:
            continue
        verts = np.stack(
            [(corners[idx] if kind == 0 else edges[idx])[sel] for kind, idx in recipe],
            axis=1,
        )
        acc += _shoelace_moments_batch(verts)
    return acc


def _leaf_corner_coords(bounds: np.ndarray) -> list[np.ndarray]:
    """The 4 corner-coordinate arrays (each ``(n, 2)``) in QuadCell.corners() order."""
    x0, x1, y0, y1 = bounds[:, 0], bounds[:, 1], bounds[:, 2], bounds[:, 3]
    return [
        np.stack([x0, y0], axis=1),
        np.stack([x1, y0], axis=1),
        np.stack([x1, y1], axis=1),
        np.stack([x0, y1], axis=1),
    ]


def _leaf_edge_crossings(corners: list[np.ndarray], phi: np.ndarray) -> list[np.ndarray]:
    """Zero-crossing point on each of the 4 edges (each ``(n, 2)``), matching
    :func:`_interpolate_zero_crossing` (degenerate edges fall back to the midpoint)."""
    edges: list[np.ndarray] = []
    for i in range(4):
        j = (i + 1) % 4
        va, vb = phi[:, i], phi[:, j]
        denom = vb - va
        degenerate = np.abs(denom) < 1e-15
        t = np.where(degenerate, 0.5, -va / np.where(degenerate, 1.0, denom))
        edges.append(corners[i] + t[:, None] * (corners[j] - corners[i]))
    return edges


