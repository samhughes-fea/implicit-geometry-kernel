"""Object-model refinement of one quadtree against several fields."""

from __future__ import annotations

import warnings

from typing import Any, Callable

import numpy as np

from implicit_geometry_kernel.integration.quadtree.refine import _active_support_mask, _cell_cached_corners, _cells_need_refine_batch, _subdivide_refined_reuse
from implicit_geometry_kernel.integration.quadtree.tree import QuadCell, QuadTreeMesh, SupportBox, _NAN4, _cell_diagonal, _cell_interior_points, _exterior_sentinel, _invalidate_leaf_cache, _n_interior_probes


def refine_quadtree_multi(
    mesh: QuadTreeMesh,
    samplers: dict[str, Any],
    supports: dict[str, SupportBox] | None = None,
    *,
    multi_sampler: Callable[[np.ndarray], dict[str, np.ndarray]] | None = None,
) -> QuadTreeMesh:
    """Refine near the union of all scalar-field zero levels.

    A cell subdivides when **any** active sampler reports a sign change or passes
    within a cell-size-relative band of zero across corner + interior probe points.
    """
    if not samplers:
        raise ValueError("samplers must be non-empty.")
    _invalidate_leaf_cache(mesh)
    spec = mesh.spec
    labels = [str(lbl) for lbl in samplers.keys()]
    sentinel = _exterior_sentinel(spec)
    n_probes = _n_interior_probes(spec)
    candidates = mesh.leaf_cells()
    leaf_count = len(candidates)
    while True:
        if leaf_count >= mesh.spec.max_cells:
            warnings.warn(
                f"refine_quadtree_multi: max_cells={mesh.spec.max_cells} reached with "
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
        n_el = len(eligible)
        cell_diags = np.array([_cell_diagonal(c) for c in eligible], dtype=float)
        active_mask = _active_support_mask(eligible, labels, supports)
        all_interior = np.vstack([_cell_interior_points(c, spec) for c in eligible])

        needs = np.zeros(n_el, dtype=bool)
        corners_by_label: dict[str, np.ndarray] = {}
        interior_by_label: dict[str, np.ndarray] = {}
        if multi_sampler is not None:
            # One batched interior-probe evaluation for all labels. The probes
            # coincide with the children's new corners, so they are reused for
            # subdivision below (see _subdivide_refined_reuse). Parent corners
            # come from the per-cell cache -- every refined-mesh cell stores
            # phi_fields -- so no separate corner evaluation is needed.
            phi_interior_batch = multi_sampler(all_interior)
            interior_true = {
                lbl: np.asarray(phi_interior_batch[lbl], dtype=float).reshape(n_el, n_probes)
                for lbl in labels
            }
            corner_cache = {
                lbl: np.full((n_el, 4), sentinel, dtype=float) for lbl in labels
            }
            missing_rows: list[int] = []
            for i, cell in enumerate(eligible):
                row_missing = False
                for j, lbl in enumerate(labels):
                    if not active_mask[i, j]:
                        continue
                    cached = _cell_cached_corners(cell, lbl, labels)
                    if cached is None:
                        row_missing = True
                    else:
                        corner_cache[lbl][i] = cached
                if row_missing:
                    missing_rows.append(i)
            if missing_rows:
                miss = np.asarray(missing_rows, dtype=int)
                miss_pts = np.vstack([eligible[i].corners() for i in miss])
                miss_batch = multi_sampler(miss_pts)
                for j, lbl in enumerate(labels):
                    vals = np.asarray(miss_batch[lbl], dtype=float).reshape(miss.size, 4)
                    active_here = active_mask[miss, j]
                    rows = miss[active_here]
                    if rows.size:
                        corner_cache[lbl][rows] = vals[active_here]
            for j, lbl in enumerate(labels):
                active_col = active_mask[:, j]
                corners = corner_cache[lbl]
                interior = np.full((n_el, n_probes), sentinel, dtype=float)
                interior[active_col] = interior_true[lbl][active_col]
                corners_by_label[lbl] = corners
                interior_by_label[lbl] = interior
                all_phi = np.concatenate([corners, interior], axis=1)
                label_needs = _cells_need_refine_batch(
                    all_phi,
                    spec.refine_band,
                    spec.refine_band_scale,
                    cell_diags,
                )
                needs |= label_needs & active_col
                sign_change = (all_phi.min(axis=1) <= 0.0) & (all_phi.max(axis=1) >= 0.0)
                mesh.refine_stats["sign_change"] += int(
                    np.sum(sign_change & label_needs & active_col)
                )
        else:
            for j, lbl in enumerate(labels):
                sampler = samplers[lbl]
                corners = np.full((n_el, 4), sentinel, dtype=float)
                interior = np.full((n_el, n_probes), sentinel, dtype=float)

                for i, cell in enumerate(eligible):
                    if not active_mask[i, j]:
                        continue
                    phi_corners = _cell_cached_corners(cell, lbl, labels)
                    if phi_corners is not None:
                        corners[i] = phi_corners

                active_rows = np.flatnonzero(active_mask[:, j])
                missing_corners = active_rows[
                    np.all(corners[active_rows] == sentinel, axis=1)
                ]
                if len(missing_corners):
                    corner_pts = np.vstack([eligible[i].corners() for i in missing_corners])
                    corner_vals = np.asarray(sampler(corner_pts), dtype=float).reshape(
                        len(missing_corners), 4
                    )
                    corners[missing_corners] = corner_vals

                if len(active_rows):
                    gather: list[int] = []
                    for i in active_rows:  # type: ignore
                        gather.extend(range(int(i) * n_probes, (int(i) + 1) * n_probes))
                    phi_int = np.asarray(sampler(all_interior[gather]), dtype=float).reshape(
                        len(active_rows), n_probes
                    )
                    interior[active_rows] = phi_int

                corners_by_label[lbl] = corners
                interior_by_label[lbl] = interior
                all_phi = np.concatenate([corners, interior], axis=1)
                label_needs = _cells_need_refine_batch(
                    all_phi,
                    spec.refine_band,
                    spec.refine_band_scale,
                    cell_diags,
                )
                needs |= label_needs & active_mask[:, j]

                sign_change = (all_phi.min(axis=1) <= 0.0) & (all_phi.max(axis=1) >= 0.0)
                mesh.refine_stats["sign_change"] += int(
                    np.sum(sign_change & label_needs & active_mask[:, j])
                )

        # Early-stop leaves: stash probe[0] (cell centre) for _attach_center_phi reuse.
        notneed_idx = np.flatnonzero(~needs)
        for i in notneed_idx:  # type: ignore
            cell = eligible[int(i)]
            cell._center_phi = {  # type: ignore[attr-defined]
                lbl: float(interior_by_label[lbl][int(i), 0]) for lbl in labels
            }

        ref_idx = np.flatnonzero(needs)
        if ref_idx.size == 0:
            break
        new_children = _subdivide_refined_reuse(
            eligible,
            ref_idx,
            labels,
            corners_by_label,
            interior_by_label,
            supports,
            spec,
        )
        mesh.refine_stats["subdivisions"] += int(ref_idx.size)
        leaf_count += 3 * int(ref_idx.size)
        candidates = new_children
    # The loop above mutated the tree (subdivisions) without touching the leaf
    # cache, so drop it here: the next leaf_cells() rebuilds from the refined
    # tree. Without this, callers that don't explicitly refresh (e.g.
    # build_area_mesh_2d) see the stale pre-refinement roots.
    _invalidate_leaf_cache(mesh)
    leaves = mesh.leaf_cells()
    center_phi: dict[str, np.ndarray] = {}
    for lbl in labels:
        vals = np.empty(len(leaves), dtype=float)  # type: ignore
        for i, cell in enumerate(leaves):
            cp = getattr(cell, "_center_phi", None)
            vals[i] = float(cp[lbl]) if cp is not None and lbl in cp else np.nan
        center_phi[lbl] = vals
    mesh._leaf_center_phi = center_phi  # type: ignore[attr-defined]
    return mesh


def _phi_by_label_from_leaves(
    leaves: list[QuadCell],
    samplers: dict[str, Any],
) -> dict[str, np.ndarray]:
    if len(samplers) == 1:
        label = next(iter(samplers.keys()))
        return {str(label): np.array([c.phi_corners for c in leaves], dtype=float)}
    return {
        lbl: np.array(
            [(c.phi_fields or {}).get(lbl, _NAN4) for c in leaves],
            dtype=float,
        )
        for lbl in samplers
    }


