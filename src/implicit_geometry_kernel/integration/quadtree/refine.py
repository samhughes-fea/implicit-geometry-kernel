"""Object-model quadtree refinement for one field, plus the multi-field cell subdivision helpers."""

from __future__ import annotations

import warnings

from typing import Any, Callable

import numpy as np

from implicit_geometry_kernel.integration.spec import QuadtreeSpec
from implicit_geometry_kernel.integration.quadtree.tree import QuadCell, QuadTreeMesh, SupportBox, _cell_diagonal, _cell_interior_points, _cell_overlaps_support, _cells_bounds_array, _exterior_sentinel, _invalidate_leaf_cache, _n_interior_probes, _sample_cell_phi, _subdivide_cell


def _cell_needs_refine(
    all_phi_cell: np.ndarray,
    band_ceiling: float,
    band_scale: float,
    cell_diag: float,
) -> bool:
    if np.min(all_phi_cell) <= 0.0 <= np.max(all_phi_cell):
        return True
    band = min(float(band_ceiling), float(band_scale) * float(cell_diag))
    return bool(np.min(np.abs(all_phi_cell)) < band)


def _cells_need_refine_batch(
    all_phi: np.ndarray,
    band_ceiling: float,
    band_scale: float,
    cell_diags: np.ndarray,
) -> np.ndarray:
    """Vectorized refine test; ``all_phi`` shape ``(n_cells, n_probes)``."""
    mins = all_phi.min(axis=1)
    maxs = all_phi.max(axis=1)
    sign_change = (mins <= 0.0) & (maxs >= 0.0)
    abs_mins = np.abs(all_phi).min(axis=1)
    bands = np.minimum(float(band_ceiling), float(band_scale) * cell_diags)
    return sign_change | (abs_mins < bands)


def _active_support_mask(
    cells: list[QuadCell],
    labels: list[str],
    supports: dict[str, SupportBox] | None,
) -> np.ndarray:
    n = len(cells)
    m = len(labels)
    if supports is None:
        return np.ones((n, m), dtype=bool)
    if n == 0:
        return np.empty((0, m), dtype=bool)
    bounds = _cells_bounds_array(cells)
    x_min = bounds[:, 0]
    x_max = bounds[:, 1]
    y_min = bounds[:, 2]
    y_max = bounds[:, 3]
    mask = np.zeros((n, m), dtype=bool)
    for j, lbl in enumerate(labels):
        box = supports.get(lbl)
        if box is None:
            mask[:, j] = True
        else:
            bx0, bx1, by0, by1 = box
            mask[:, j] = ~(
                (x_max <= bx0)
                | (x_min >= bx1)
                | (y_max <= by0)
                | (y_min >= by1)
            )
    return mask


def refine_quadtree(mesh: QuadTreeMesh, phi_sampler) -> QuadTreeMesh:
    """Iteratively refine the mesh near the scalar-field zero level."""
    _invalidate_leaf_cache(mesh)
    spec = mesh.spec
    n_probes = _n_interior_probes(spec)
    candidates = mesh.leaf_cells()
    leaf_count = len(candidates)
    while True:
        if leaf_count >= mesh.spec.max_cells:
            warnings.warn(
                f"refine_quadtree: max_cells={mesh.spec.max_cells} reached with "
                f"{leaf_count} leaves; mesh is truncated.",
                stacklevel=2,
            )
            break
        eligible = [
            c
            for c in candidates
            if c.level < mesh.spec.max_depth
            and max(c.width, c.height) > mesh.spec.min_cell_size
        ]
        if not eligible:
            break
        all_interior = np.vstack([_cell_interior_points(c, spec) for c in eligible])
        all_phi = np.asarray(phi_sampler(all_interior), dtype=float)
        cell_diags = np.array([_cell_diagonal(c) for c in eligible], dtype=float)
        to_refine = []
        for i, cell in enumerate(eligible):
            phi_corners = (
                cell.phi_corners
                if cell.phi_corners is not None
                else _sample_cell_phi(cell, phi_sampler)
            )
            phi_interior = all_phi[i * n_probes : (i + 1) * n_probes]
            all_phi_cell = np.concatenate([phi_corners, phi_interior])
            if _cell_needs_refine(
                all_phi_cell,
                spec.refine_band,
                spec.refine_band_scale,
                cell_diags[i],
            ):
                to_refine.append(cell)
                if np.min(all_phi_cell) <= 0.0 <= np.max(all_phi_cell):
                    mesh.refine_stats["sign_change"] += 1
                else:
                    mesh.refine_stats["near_interface"] += 1
        if not to_refine:
            break
        new_children: list[QuadCell] = []
        for cell in to_refine:
            children = _subdivide_cell(cell, phi_sampler)
            new_children.extend(children)
            mesh.refine_stats["subdivisions"] += 1
        leaf_count += 3 * len(to_refine)
        candidates = new_children
    # The loop above mutated the tree (subdivisions) without touching the leaf
    # cache, so drop it here: the next leaf_cells() rebuilds from the refined
    # tree. Without this, callers that don't explicitly refresh (e.g.
    # build_area_mesh_2d) see the stale pre-refinement roots.
    _invalidate_leaf_cache(mesh)
    return mesh


def _cell_cached_corners(
    cell: QuadCell,
    lbl: str,
    labels: list[str],
) -> np.ndarray | None:
    if cell.phi_fields is not None and lbl in cell.phi_fields:
        return cell.phi_fields[lbl]
    if cell.phi_corners is not None and len(labels) == 1:
        return cell.phi_corners
    return None


def _subdivide_cell_multi(
    cell: QuadCell,
    samplers: dict[str, Any],
    *,
    supports: dict[str, SupportBox] | None = None,
    spec: QuadtreeSpec | None = None,
    multi_sampler: Callable[[np.ndarray], dict[str, np.ndarray]] | None = None,
) -> list[QuadCell]:
    """Subdivide a cell and sample active field samplers at child corners."""
    if spec is None:
        raise ValueError("spec is required for multi-field subdivision.")
    sentinel = _exterior_sentinel(spec)
    xm = 0.5 * (cell.x_min + cell.x_max)
    ym = 0.5 * (cell.y_min + cell.y_max)
    children = [
        QuadCell(cell.x_min, xm, cell.y_min, ym, level=cell.level + 1),
        QuadCell(xm, cell.x_max, cell.y_min, ym, level=cell.level + 1),
        QuadCell(xm, cell.x_max, ym, cell.y_max, level=cell.level + 1),
        QuadCell(cell.x_min, xm, ym, cell.y_max, level=cell.level + 1),
    ]
    all_corners = np.vstack([c.corners() for c in children])
    phi_fields: dict[str, np.ndarray] = {}
    if multi_sampler is not None:
        batch = multi_sampler(all_corners)
        for label in samplers:
            lbl = str(label)
            vals = np.asarray(batch[lbl], dtype=float).reshape(-1)
            if supports is not None and lbl in supports:
                support = supports[lbl]
                for i, child in enumerate(children):
                    if not _cell_overlaps_support(child, support):
                        vals[i * 4 : (i + 1) * 4] = sentinel
            phi_fields[lbl] = vals
    else:
        for label, sampler in samplers.items():
            lbl = str(label)
            vals = np.asarray(sampler(all_corners), dtype=float)
            if supports is not None and lbl in supports:
                support = supports[lbl]
                for i, child in enumerate(children):
                    if not _cell_overlaps_support(child, support):
                        vals[i * 4 : (i + 1) * 4] = sentinel
            phi_fields[lbl] = vals
    for i, child in enumerate(children):
        child.phi_fields = {
            lbl: vals[i * 4 : (i + 1) * 4] for lbl, vals in phi_fields.items()
        }
        if len(phi_fields) == 1:
            child.phi_corners = next(iter(child.phi_fields.values()))
    cell.children = children
    return children


# Maps each child's 4 corners (in QuadCell.corners() order) to a source in the
# parent's decision-phase samples: ("c", k) = parent corner k, ("p", k) =
# interior probe k. The default 5 probes from _cell_interior_points are
# [center, bottom-mid, right-mid, top-mid, left-mid] -- exactly the new vertices
# subdivision introduces -- so child-corner phi is a gather over already-sampled
# values rather than a fresh evaluation.
_CHILD_CORNER_RECIPE: tuple[tuple[tuple[str, int], ...], ...] = (
    (("c", 0), ("p", 1), ("p", 0), ("p", 4)),  # child 0: (x_min, xm, y_min, ym)
    (("p", 1), ("c", 1), ("p", 2), ("p", 0)),  # child 1: (xm, x_max, y_min, ym)
    (("p", 0), ("p", 2), ("c", 2), ("p", 3)),  # child 2: (xm, x_max, ym, y_max)
    (("p", 4), ("p", 0), ("p", 3), ("c", 3)),  # child 3: (x_min, xm, ym, y_max)
)


def _subdivide_refined_reuse(
    eligible: list[QuadCell],
    ref_idx: np.ndarray,
    labels: list[str],
    corners_by_label: dict[str, np.ndarray],
    interior_by_label: dict[str, np.ndarray],
    supports: dict[str, SupportBox] | None,
    spec: QuadtreeSpec,
) -> list[QuadCell]:
    """Subdivide refined cells, reusing decision-phase corner + probe phi.

    The refine decision already evaluated every parent's 4 corners and 5 interior
    probes, which coincide with the corners of the 4 children (see
    ``_CHILD_CORNER_RECIPE``). Child-corner phi is therefore a gather over those
    arrays -- no new field evaluation -- and per-child support sentinel masking
    reproduces :func:`_subdivide_cell_multi` bit-for-bit.
    """
    sentinel = _exterior_sentinel(spec)
    single = len(labels) == 1
    ref = np.asarray(ref_idx, dtype=int)
    r_count = int(ref.size)

    px0 = np.array([eligible[i].x_min for i in ref], dtype=float)
    px1 = np.array([eligible[i].x_max for i in ref], dtype=float)
    py0 = np.array([eligible[i].y_min for i in ref], dtype=float)
    py1 = np.array([eligible[i].y_max for i in ref], dtype=float)
    pxm = 0.5 * (px0 + px1)
    pym = 0.5 * (py0 + py1)
    # Child bounds (x_min, x_max, y_min, y_max), each component shape (r_count,).
    child_bounds = (
        (px0, pxm, py0, pym),  # child 0
        (pxm, px1, py0, pym),  # child 1
        (pxm, px1, pym, py1),  # child 2
        (px0, pxm, pym, py1),  # child 3
    )

    # Per label, build the 4 children's corner phi arrays (each (r_count, 4)).
    child_phi: dict[str, list[np.ndarray]] = {}
    for lbl in labels:
        c = corners_by_label[lbl][ref]      # (r_count, 4)
        p = interior_by_label[lbl][ref]     # (r_count, n_probes)
        box = supports.get(lbl) if supports is not None else None
        per_child: list[np.ndarray] = []
        for k in range(4):
            cols = [
                (c[:, idx] if kind == "c" else p[:, idx])
                for kind, idx in _CHILD_CORNER_RECIPE[k]
            ]
            vals = np.stack(cols, axis=1)   # (r_count, 4)
            if box is not None:
                bx0, bx1, by0, by1 = box
                cx0, cx1, cy0, cy1 = child_bounds[k]
                overlap = ~(
                    (cx1 <= bx0) | (cx0 >= bx1) | (cy1 <= by0) | (cy0 >= by1)
                )
                vals[~overlap] = sentinel
            per_child.append(vals)
        child_phi[lbl] = per_child

    new_children: list[QuadCell] = []
    for r in range(r_count):
        parent = eligible[int(ref[r])]
        lvl = parent.level + 1
        children = [
            QuadCell(
                float(child_bounds[k][0][r]),
                float(child_bounds[k][1][r]),
                float(child_bounds[k][2][r]),
                float(child_bounds[k][3][r]),
                level=lvl,
            )
            for k in range(4)
        ]
        for k, child in enumerate(children):
            fields = {lbl: child_phi[lbl][k][r] for lbl in labels}
            child.phi_fields = fields
            if single:
                child.phi_corners = next(iter(fields.values()))
        parent.children = children
        new_children.extend(children)
    return new_children


