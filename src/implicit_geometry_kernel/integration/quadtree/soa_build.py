"""SoA refinement loop and shared multi-field quadtree."""

from __future__ import annotations

import warnings
from dataclasses import replace
from typing import Any, Callable

import numpy as np

from implicit_geometry_kernel.integration.spec import QuadtreeSpec
from implicit_geometry_kernel.integration.quadtree.tree import (
    QuadCell,
    QuadTreeMesh,
    SupportBox,
    _exterior_sentinel,
    _n_interior_probes,
    _square_forest_domain,
)
from implicit_geometry_kernel.integration.quadtree.soa import (
    BatchEval,
    _active_support_mask_bounds,
    _cells_need_refine_batch,
    _expand_children_soa,
    _refine_points_for_bounds,
    _samplers_as_batch,
)


def _refine_soa(
    spec: QuadtreeSpec,
    labels: list[str],
    batch_eval: BatchEval,
    supports: dict[str, SupportBox] | None,
    frontier: tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]],
) -> QuadTreeMesh:
    """Struct-of-arrays adaptive refinement (no per-cell ``QuadCell`` objects).

    Mirrors the SoA refine decision (corner cache + interior probe
    band test, per-label support culling) but carries the frontier as arrays and
    expands children by reshaping. Produces the same leaf *set* and count as the
    object path; leaf *order* is refinement/level order (summation reorders, so
    properties match to float round-off — see plan parity note).
    """
    sentinel = _exterior_sentinel(spec)
    n_probes = _n_interior_probes(spec)
    max_depth = spec.max_depth
    min_cell = spec.min_cell_size
    max_cells = spec.max_cells
    band_ceiling = spec.refine_band
    band_scale = spec.refine_band_scale
    thin_band = bool(spec.thin_band)

    fb, flevel, fcorners = frontier
    refine_stats = {"sign_change": 0, "near_interface": 0, "subdivisions": 0}

    lb_chunks: list[np.ndarray] = []
    lv_chunks: list[np.ndarray] = []
    lc_chunks: dict[str, list[np.ndarray]] = {lbl: [] for lbl in labels}
    # Centre φ stash: finite for early-stop leaves (probe[0]); NaN = must eval in attach.
    lctr_chunks: dict[str, list[np.ndarray]] = {lbl: [] for lbl in labels}

    def _emit(mask, b, lev, cor, centers=None):
        lb_chunks.append(b[mask])
        lv_chunks.append(lev[mask])
        k = int(np.count_nonzero(mask))
        for lbl in labels:
            lc_chunks[lbl].append(cor[lbl][mask])
            if centers is None:
                lctr_chunks[lbl].append(np.full(k, np.nan, dtype=float))
            else:
                lctr_chunks[lbl].append(np.asarray(centers[lbl][mask], dtype=float))

    leaf_count = fb.shape[0]
    while True:
        if fb.shape[0] == 0:
            break
        if leaf_count >= max_cells:
            warnings.warn(
                f"_refine_soa: max_cells={max_cells} reached with {leaf_count} "
                f"leaves; mesh is truncated.",
                stacklevel=2,
            )
            _emit(np.ones(fb.shape[0], dtype=bool), fb, flevel, fcorners)
            break

        w = fb[:, 1] - fb[:, 0]
        h = fb[:, 3] - fb[:, 2]
        eligible = (flevel < max_depth) & (np.maximum(w, h) > min_cell)
        if not eligible.any():
            _emit(np.ones(fb.shape[0], dtype=bool), fb, flevel, fcorners)
            break

        inel = ~eligible
        if inel.any():
            _emit(inel, fb, flevel, fcorners)

        eb = fb[eligible]
        elev = flevel[eligible]
        ecorners = {lbl: fcorners[lbl][eligible] for lbl in labels}
        n_el = eb.shape[0]

        probe_pts = _refine_points_for_bounds(eb, spec).reshape(-1, 2)
        batch = batch_eval(probe_pts)
        interior_true = {
            lbl: np.asarray(batch[lbl], dtype=float).reshape(n_el, n_probes) for lbl in labels
        }
        cell_diags = np.hypot(eb[:, 1] - eb[:, 0], eb[:, 3] - eb[:, 2])
        active = _active_support_mask_bounds(eb, labels, supports)

        needs = np.zeros(n_el, dtype=bool)
        for j, lbl in enumerate(labels):
            acol = active[:, j]
            interior = interior_true[lbl].copy()
            interior[~acol] = sentinel
            all_phi = np.concatenate([ecorners[lbl], interior], axis=1)
            label_needs = _cells_need_refine_batch(
                all_phi, band_ceiling, band_scale, cell_diags, thin_band
            )
            needs |= label_needs & acol
            sign_change = (all_phi.min(axis=1) <= 0.0) & (all_phi.max(axis=1) >= 0.0)
            refine_stats["sign_change"] += int(np.sum(sign_change & label_needs & acol))

        notneed = ~needs
        if notneed.any():
            # probe[0] is the cell centre — reuse in _attach_center_phi.
            centers = {}
            for j, lbl in enumerate(labels):
                cvals = np.asarray(interior_true[lbl][:, 0], dtype=float).copy()
                cvals[~active[:, j]] = sentinel
                centers[lbl] = cvals
            _emit(notneed, eb, elev, ecorners, centers)

        ref = needs
        R = int(np.count_nonzero(ref))
        if R == 0:
            break

        rb = eb[ref]
        rlev = elev[ref]
        rcorners = {lbl: ecorners[lbl][ref] for lbl in labels}
        rprobes = {lbl: interior_true[lbl][ref] for lbl in labels}
        fb, flevel, fcorners = _expand_children_soa(
            rb, rlev, rcorners, rprobes, labels, supports, spec
        )
        refine_stats["subdivisions"] += R
        leaf_count += 3 * R

    if lb_chunks:
        leaf_bounds_arr = np.concatenate(lb_chunks, axis=0)
        leaf_levels_arr = np.concatenate(lv_chunks, axis=0).astype(int)
        phi_by_label = {lbl: np.concatenate(lc_chunks[lbl], axis=0) for lbl in labels}
        leaf_center_phi = {lbl: np.concatenate(lctr_chunks[lbl], axis=0) for lbl in labels}
    else:
        leaf_bounds_arr = np.empty((0, 4), dtype=float)
        leaf_levels_arr = np.empty((0,), dtype=int)
        phi_by_label = {lbl: np.empty((0, 4), dtype=float) for lbl in labels}
        leaf_center_phi = {lbl: np.empty((0,), dtype=float) for lbl in labels}

    mesh = QuadTreeMesh(spec=spec, roots=[])
    mesh.refine_stats = refine_stats
    mesh._is_soa = True  # type: ignore[attr-defined]
    mesh._leaf_bounds = leaf_bounds_arr  # type: ignore[attr-defined]
    mesh._leaf_levels = leaf_levels_arr  # type: ignore[attr-defined]
    mesh._cached_leaves = None  # type: ignore[attr-defined]
    mesh.phi_by_label = phi_by_label  # type: ignore[attr-defined]
    mesh._leaf_center_phi = leaf_center_phi  # type: ignore[attr-defined]
    return mesh


def _soa_frontier_corners(
    bounds: np.ndarray, labels: list[str], batch_eval: BatchEval
) -> dict[str, np.ndarray]:
    """Evaluate root-cell corner φ for the initial SoA frontier."""
    m = bounds.shape[0]
    x0, x1, y0, y1 = bounds[:, 0], bounds[:, 1], bounds[:, 2], bounds[:, 3]
    # Corner order matches QuadCell.corners(): (x0,y0),(x1,y0),(x1,y1),(x0,y1).
    corner_pts = np.empty((m, 4, 2), dtype=float)
    corner_pts[:, 0] = np.stack([x0, y0], axis=1)
    corner_pts[:, 1] = np.stack([x1, y0], axis=1)
    corner_pts[:, 2] = np.stack([x1, y1], axis=1)
    corner_pts[:, 3] = np.stack([x0, y1], axis=1)
    batch = batch_eval(corner_pts.reshape(-1, 2))
    return {lbl: np.asarray(batch[lbl], dtype=float).reshape(m, 4) for lbl in labels}


def _square_forest_roots(spec: QuadtreeSpec) -> tuple[QuadtreeSpec, np.ndarray]:
    """Square-tiled root-cell bounds ``(n_x*n_y, 4)`` and the expanded forest spec."""
    x_min, x_max, y_min, y_max, s, n_x, n_y = _square_forest_domain(spec)
    forest_spec = replace(
        spec, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, root_cell_size=s
    )
    rows = []
    for j in range(n_y):
        for i in range(n_x):
            rows.append([x_min + i * s, x_min + (i + 1) * s, y_min + j * s, y_min + (j + 1) * s])
    return forest_spec, np.array(rows, dtype=float)


def build_shared_quadtree(
    spec: QuadtreeSpec,
    samplers: dict[str, Any],
    supports: dict[str, SupportBox] | None = None,
    *,
    multi_sampler: Callable[[np.ndarray], dict[str, np.ndarray]] | None = None,
) -> QuadTreeMesh:
    """Build and refine one quadtree (or square-root forest) for named scalar fields.

    Struct-of-arrays refinement (``_refine_soa``): the frontier is carried as
    arrays and children are expanded by reshaping, with no per-cell ``QuadCell``
    objects in the loop. The returned mesh exposes ``leaf_bounds`` /
    ``phi_by_label`` directly; ``leaf_cells()`` materialises objects lazily for
    consumers that still iterate cells.
    """
    if not samplers:
        raise ValueError("samplers must be non-empty.")
    labels = [str(lbl) for lbl in samplers.keys()]
    batch_eval = multi_sampler if multi_sampler is not None else _samplers_as_batch(samplers, labels)

    if spec.square_cells:
        run_spec, root_bounds = _square_forest_roots(spec)
    else:
        run_spec = spec
        root_bounds = np.array(
            [[spec.x_min, spec.x_max, spec.y_min, spec.y_max]], dtype=float
        )

    root_corners = _soa_frontier_corners(root_bounds, labels, batch_eval)
    root_levels = np.zeros(root_bounds.shape[0], dtype=int)
    mesh = _refine_soa(run_spec, labels, batch_eval, supports, (root_bounds, root_levels, root_corners))

    # Keep the root cells for introspection (e.g. square-vs-rectangular checks);
    # leaves come from the SoA arrays, not this (childless) root list.
    mesh.roots = [
        QuadCell(float(b[0]), float(b[1]), float(b[2]), float(b[3]), level=0)
        for b in root_bounds
    ]
    if supports is not None:
        mesh.component_supports = supports  # type: ignore[attr-defined]
    return mesh


