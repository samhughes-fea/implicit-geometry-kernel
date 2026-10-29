"""Section property quadrature on adaptive quadtree leaf cells."""

from __future__ import annotations

from typing import Any, Callable

import numpy as np

from .quadtree import (
    QuadCell,
    QuadTreeMesh,
    SupportBox,
    _exterior_sentinel,
    _interpolate_zero_crossing,
    leaf_bounds,
)


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
    edge_pairs = [(0, 1), (1, 2), (2, 3), (3, 0)]
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


def _label_corner_phi(mesh: QuadTreeMesh, label: Any) -> np.ndarray:
    key = str(getattr(label, "value", label))
    return np.asarray(mesh.phi_by_label[key], dtype=float)[:, :4]


def label_fill_polygons(mesh: QuadTreeMesh, label: Any) -> list[np.ndarray]:
    """Filled ``φ <= 0`` polygons (one per non-exterior leaf) for ``label``.

    Interior leaves (all corners ≤ 0) contribute their full rectangle; boundary
    leaves contribute the marching-squares inside-polygon. The union tiles the
    ``φ < 0`` region with no gaps, so it is a watertight fill for rendering. Uses
    the same per-cell edge-crossing recipe as the property quadrature, so the
    summed polygon area matches the integrated section area.
    """
    phi = _label_corner_phi(mesh, label)
    bounds = leaf_bounds(mesh)
    n = bounds.shape[0]
    polys: list[np.ndarray] = []
    if n == 0:
        return polys
    corners = _leaf_corner_coords(bounds)
    inside = phi <= 0.0
    codes = (
        inside[:, 0].astype(np.int64)
        | (inside[:, 1].astype(np.int64) << 1)
        | (inside[:, 2].astype(np.int64) << 2)
        | (inside[:, 3].astype(np.int64) << 3)
    )
    edges = _leaf_edge_crossings(corners, phi)
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
        )  # (m, L, 2)
        polys.extend(verts)
    return polys


def label_boundary_segments(mesh: QuadTreeMesh, label: Any) -> np.ndarray:
    """φ=0 marching-squares segments ``(m, 2, 2)`` for ``label`` (unstitched).

    One segment per 2-crossing boundary leaf, two per 4-crossing (saddle) leaf,
    matching :func:`contour_from_quadtree_marching_squares`. Suitable for a
    ``LineCollection`` boundary overlay; minor gaps at adaptive T-junctions are
    cosmetically negligible at plot scale.
    """
    phi = _label_corner_phi(mesh, label)
    bounds = leaf_bounds(mesh)
    n = bounds.shape[0]
    if n == 0:
        return np.empty((0, 2, 2), dtype=float)
    inside = phi <= 0.0
    n_in = inside.sum(axis=1)
    boundary = np.flatnonzero((n_in > 0) & (n_in < 4))
    if boundary.size == 0:
        return np.empty((0, 2, 2), dtype=float)
    corners = _leaf_corner_coords(bounds)
    edges = _leaf_edge_crossings(corners, phi)
    crossing = np.column_stack([inside[:, i] != inside[:, (i + 1) % 4] for i in range(4)])
    segs: list[list[np.ndarray]] = []
    for c in boundary:
        pts = [edges[e][c] for e in range(4) if crossing[c, e]]
        if len(pts) == 2:
            segs.append([pts[0], pts[1]])
        elif len(pts) == 4:
            segs.append([pts[0], pts[1]])
            segs.append([pts[2], pts[3]])
    return np.asarray(segs, dtype=float) if segs else np.empty((0, 2, 2), dtype=float)


def _accumulate_moments_vectorized(
    acc: list[float],
    bounds: np.ndarray,
    phi_samples: np.ndarray,
    *,
    sampler: Callable[[np.ndarray], np.ndarray] | None = None,
    max_subdiv: int = 2,
) -> None:
    """Classify leaves in batch and accumulate interior/boundary/hollow contributions.

    Operates purely on the ``bounds`` ``(n, 4)`` and per-leaf φ arrays, so it
    needs no materialised ``QuadCell`` list (the rare hollow branch builds a
    temporary cell from the leaf's bounds).
    """
    n = bounds.shape[0]
    if n == 0:
        return
    phi = np.asarray(phi_samples, dtype=float).reshape(n, -1)
    corners = phi[:, :4]
    if phi.shape[1] >= 5:
        center = phi[:, 4]
    else:
        center = corners.mean(axis=1)

    exterior = np.all(corners > 0.0, axis=1) & (center > 0.0)
    interior = np.all(corners <= 0.0, axis=1) & (center <= 0.0)
    hollow = np.all(corners > 0.0, axis=1) & (center <= 0.0)
    boundary = ~(exterior | interior | hollow)

    if np.any(interior):
        idx = np.flatnonzero(interior)
        b = bounds[idx]
        moments = _rectangle_moments_batch(b[:, 0], b[:, 1], b[:, 2], b[:, 3])
        acc_arr = np.sum(moments, axis=0)
        for i in range(6):
            acc[i] += float(acc_arr[i])

    if np.any(boundary):
        idx = np.flatnonzero(boundary)
        acc_b = _boundary_moments_batch(bounds[idx], corners[idx])
        for j in range(6):
            acc[j] += float(acc_b[j])

    if np.any(hollow) and sampler is not None:
        for i in np.flatnonzero(hollow):
            b = bounds[i]
            cell = QuadCell(float(b[0]), float(b[1]), float(b[2]), float(b[3]))
            _accumulate_cell(
                acc,
                cell,
                phi[i],
                sampler=sampler,
                subdiv_depth=0,
                max_subdiv=max_subdiv,
            )


def _raw_to_section_properties(acc: list[float]) -> dict[str, float]:
    area, sx, sy, sxx, syy, sxy = acc
    if area <= 0.0:
        nan = float("nan")
        return {
            "area": 0.0,
            "cx": nan,
            "cy": nan,
            "Ixx": nan,
            "Iyy": nan,
            "Ixy": nan,
            "r_gyr_x": nan,
            "r_gyr_y": nan,
        }
    cx = sx / area
    cy = sy / area
    ixx = syy - sy * sy / area
    iyy = sxx - sx * sx / area
    ixy = sxy - sx * sy / area
    return {
        "area": float(area),
        "cx": float(cx),
        "cy": float(cy),
        "Ixx": float(ixx),
        "Iyy": float(iyy),
        "Ixy": float(ixy),
        "r_gyr_x": float(np.sqrt(ixx / area)) if area > 0 else float("nan"),
        "r_gyr_y": float(np.sqrt(iyy / area)) if area > 0 else float("nan"),
    }


def _active_leaf_indices(
    bounds: np.ndarray,
    supports: dict[str, SupportBox] | None,
    lbl: str,
) -> np.ndarray:
    """Leaf indices whose bounds overlap ``lbl`` support (or all if unbounded)."""
    n_leaves = bounds.shape[0]
    if supports is None:
        return np.arange(n_leaves, dtype=int)
    box = supports.get(lbl)
    if box is None:
        return np.arange(n_leaves, dtype=int)
    bx0, bx1, by0, by1 = box
    return np.flatnonzero(
        ~(
            (bounds[:, 1] <= bx0)
            | (bounds[:, 0] >= bx1)
            | (bounds[:, 3] <= by0)
            | (bounds[:, 2] >= by1)
        )
    )


def _attach_center_phi(
    mesh: QuadTreeMesh,
    phi_by_label: dict[str, np.ndarray],
    samplers: dict[str, Callable[[np.ndarray], np.ndarray]],
    supports: dict[str, SupportBox] | None = None,
    *,
    multi_sampler: Callable[[np.ndarray], dict[str, np.ndarray]] | None = None,
) -> dict[str, np.ndarray]:
    """Append cell-centre φ samples (shape ``(n_leaves, 5)``).

    When ``mesh._leaf_center_phi`` is present (SoA / object refine stash), finite
    entries are reused bit-identically and only NaN (missing) centres are
    evaluated — typically max-depth interface leaves that never had probe[0].
    """
    bounds = leaf_bounds(mesh)
    n_leaves = bounds.shape[0]
    if n_leaves == 0:
        return phi_by_label
    if supports is None:
        supports = getattr(mesh, "component_supports", None)
    centers = np.column_stack(
        [
            0.5 * (bounds[:, 0] + bounds[:, 1]),
            0.5 * (bounds[:, 2] + bounds[:, 3]),
        ]
    )
    sentinel = _exterior_sentinel(mesh.spec)
    stashed: dict[str, np.ndarray] | None = getattr(mesh, "_leaf_center_phi", None)

    missing_any = np.ones(n_leaves, dtype=bool)
    if stashed is not None:
        missing_any = np.zeros(n_leaves, dtype=bool)
        for lbl in samplers:
            s = stashed.get(lbl)
            if s is None or np.asarray(s).shape[0] != n_leaves:
                missing_any[:] = True
                break
            missing_any |= np.isnan(np.asarray(s, dtype=float))

    miss_idx = np.flatnonzero(missing_any)
    eval_on_miss: dict[str, np.ndarray] = {}
    if miss_idx.size:
        if multi_sampler is not None:
            batch = multi_sampler(centers[miss_idx])
            eval_on_miss = {
                lbl: np.asarray(batch[lbl], dtype=float).reshape(-1) for lbl in batch
            }
        else:
            for lbl, sampler in samplers.items():
                if sampler is None:
                    continue
                active_idx = _active_leaf_indices(bounds, supports, lbl)
                active_miss = active_idx[np.isin(active_idx, miss_idx)]
                vals = np.full(miss_idx.size, sentinel, dtype=float)
                if active_miss.size:
                    pos = np.searchsorted(miss_idx, active_miss)
                    vals[pos] = np.asarray(
                        sampler(centers[active_miss]), dtype=float
                    ).reshape(-1)
                eval_on_miss[lbl] = vals

    out: dict[str, np.ndarray] = {}
    for lbl, corners in phi_by_label.items():
        if lbl not in samplers and lbl not in eval_on_miss:
            if stashed is None or lbl not in stashed:
                out[lbl] = corners
                continue
        center_col = np.full(n_leaves, sentinel, dtype=float)
        active_idx = _active_leaf_indices(bounds, supports, lbl)
        active_mask = np.zeros(n_leaves, dtype=bool)
        active_mask[active_idx] = True

        s_lbl = None
        if stashed is not None and lbl in stashed and np.asarray(stashed[lbl]).shape[0] == n_leaves:
            s_lbl = np.asarray(stashed[lbl], dtype=float)
            known = ~np.isnan(s_lbl)
            # Stash already applied support→sentinel for early-stop leaves.
            center_col[known] = s_lbl[known]

        if lbl in eval_on_miss and miss_idx.size:
            # Only fill leaves that were missing for *this* label (another label
            # being missing can put the leaf in miss_idx without invalidating stash).
            if s_lbl is not None:
                need = np.isnan(s_lbl[miss_idx]) & active_mask[miss_idx]
            else:
                need = active_mask[miss_idx]
            center_col[miss_idx[need]] = eval_on_miss[lbl][need]

        out[lbl] = np.hstack([corners, center_col.reshape(-1, 1)])
    return out


def section_properties_from_quadtree(
    mesh: QuadTreeMesh,
    phi_by_label: dict[str, np.ndarray],
    label: str,
    *,
    sampler: Callable[[np.ndarray], np.ndarray] | None = None,
) -> dict[str, float]:
    """Integrate section properties for one component over shared quadtree leaves."""
    key = str(getattr(label, "value", label))
    phi_leaf = np.asarray(phi_by_label[key], dtype=float)
    bounds = leaf_bounds(mesh)
    acc = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    _accumulate_moments_vectorized(acc, bounds, phi_leaf, sampler=sampler)
    return _raw_to_section_properties(acc)


def section_properties_all_labels(
    mesh: QuadTreeMesh,
    phi_by_label: dict[str, np.ndarray],
    labels: list[str],
    *,
    samplers: dict[str, Callable[[np.ndarray], np.ndarray]] | None = None,
) -> dict[str, dict[str, float]]:
    """Compute properties for every label on the shared mesh."""
    samplers = samplers or {}
    return {
        label: section_properties_from_quadtree(
            mesh,
            phi_by_label,
            label,
            sampler=samplers.get(str(getattr(label, "value", label))),
        )
        for label in labels
    }
