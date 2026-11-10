"""Array-based (structure-of-arrays) quadtree refinement and the shared multi-field quadtree."""

from __future__ import annotations


from typing import Any, Callable

import numpy as np

from implicit_geometry_kernel.integration.spec import QuadtreeSpec
from implicit_geometry_kernel.integration.quadtree.tree import SupportBox, _exterior_sentinel


BatchEval = Callable[[np.ndarray], dict[str, np.ndarray]]

# Maps each child's 4 corners (in QuadCell.corners() order) to a source in the
# parent's decision-phase samples: ("c", k) = parent corner k, ("p", k) =
# interior probe k.
_CHILD_CORNER_RECIPE: tuple[tuple[tuple[str, int], ...], ...] = (
    (("c", 0), ("p", 1), ("p", 0), ("p", 4)),  # child 0: (x_min, xm, y_min, ym)
    (("p", 1), ("c", 1), ("p", 2), ("p", 0)),  # child 1: (xm, x_max, y_min, ym)
    (("p", 0), ("p", 2), ("c", 2), ("p", 3)),  # child 2: (xm, x_max, ym, y_max)
    (("p", 4), ("p", 0), ("p", 3), ("c", 3)),  # child 3: (x_min, xm, ym, y_max)
)


def _effective_refine_bands(
    band_ceiling: float,
    band_scale: float,
    cell_diags: np.ndarray,
    thin_band: bool,
) -> np.ndarray:
    """Per-cell refine radius.

    Filled solids cap the radius at ``refine_band`` so cavities stay coarse.
    Thin tubes (``thin_band``) use a Lipschitz bound ``scale * diag + radius``
    so a coarse cell whose probes sit in the cavity still sees the wall.
    """
    scaled = float(band_scale) * np.asarray(cell_diags, dtype=float)
    if thin_band:
        return scaled + float(band_ceiling)
    return np.minimum(float(band_ceiling), scaled)


def _cell_needs_refine(
    all_phi_cell: np.ndarray,
    band_ceiling: float,
    band_scale: float,
    cell_diag: float,
    thin_band: bool = False,
) -> bool:
    if np.min(all_phi_cell) <= 0.0 <= np.max(all_phi_cell):
        return True
    band = float(_effective_refine_bands(band_ceiling, band_scale, np.array([cell_diag]), thin_band)[0])
    return bool(np.min(np.abs(all_phi_cell)) < band)


def _cells_need_refine_batch(
    all_phi: np.ndarray,
    band_ceiling: float,
    band_scale: float,
    cell_diags: np.ndarray,
    thin_band: bool = False,
) -> np.ndarray:
    """Vectorized refine test; ``all_phi`` shape ``(n_cells, n_probes)``."""
    mins = all_phi.min(axis=1)
    maxs = all_phi.max(axis=1)
    sign_change = (mins <= 0.0) & (maxs >= 0.0)
    abs_mins = np.abs(all_phi).min(axis=1)
    bands = _effective_refine_bands(band_ceiling, band_scale, cell_diags, thin_band)
    return sign_change | (abs_mins < bands)


def _samplers_as_batch(samplers: dict[str, Any], labels: list[str]) -> BatchEval:
    """Wrap per-label samplers as a single batch evaluator (no shared memo)."""

    def batch(pts: np.ndarray) -> dict[str, np.ndarray]:
        return {lbl: np.asarray(samplers[lbl](pts), dtype=float).reshape(-1) for lbl in labels}

    return batch


def _refine_points_for_bounds(bounds: np.ndarray, spec: QuadtreeSpec) -> np.ndarray:
    """Interior probe points ``(M, n_probes, 2)`` for a batch of cells.

    Column order matches :func:`_cell_interior_points` exactly so that the first
    five probes (centre + 4 edge midpoints) coincide with the children's new
    corners (see ``_CHILD_CORNER_RECIPE``).
    """
    x0, x1, y0, y1 = bounds[:, 0], bounds[:, 1], bounds[:, 2], bounds[:, 3]
    xm = 0.5 * (x0 + x1)
    ym = 0.5 * (y0 + y1)
    cols = [
        np.stack([xm, ym], axis=1),
        np.stack([xm, y0], axis=1),
        np.stack([x1, ym], axis=1),
        np.stack([xm, y1], axis=1),
        np.stack([x0, ym], axis=1),
    ]
    if str(spec.refine_probes).lower() == "full":
        w = x1 - x0
        h = y1 - y0
        xq1 = x0 + 0.25 * w
        xq3 = x0 + 0.75 * w
        yq1 = y0 + 0.25 * h
        yq3 = y0 + 0.75 * h
        cols += [
            np.stack([xq1, y0], axis=1),
            np.stack([xq3, y0], axis=1),
            np.stack([x1, yq1], axis=1),
            np.stack([x1, yq3], axis=1),
            np.stack([xq1, y1], axis=1),
            np.stack([xq3, y1], axis=1),
            np.stack([x0, yq1], axis=1),
            np.stack([x0, yq3], axis=1),
        ]
    return np.stack(cols, axis=1)


def _active_support_mask_bounds(
    bounds: np.ndarray,
    labels: list[str],
    supports: dict[str, SupportBox] | None,
) -> np.ndarray:
    """Vectorized per-label support overlap mask ``(n, n_labels)`` from a bounds array."""
    n = bounds.shape[0]
    m = len(labels)
    if supports is None:
        return np.ones((n, m), dtype=bool)
    if n == 0:
        return np.empty((0, m), dtype=bool)
    x_min, x_max, y_min, y_max = bounds[:, 0], bounds[:, 1], bounds[:, 2], bounds[:, 3]
    mask = np.zeros((n, m), dtype=bool)
    for j, lbl in enumerate(labels):
        box = supports.get(lbl)
        if box is None:
            mask[:, j] = True
        else:
            bx0, bx1, by0, by1 = box
            mask[:, j] = ~(
                (x_max <= bx0) | (x_min >= bx1) | (y_max <= by0) | (y_min >= by1)
            )
    return mask


def _expand_children_soa(
    bounds: np.ndarray,
    levels: np.ndarray,
    corners: dict[str, np.ndarray],
    probes: dict[str, np.ndarray],
    labels: list[str],
    supports: dict[str, SupportBox] | None,
    spec: QuadtreeSpec,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Vectorized subdivision of ``R`` refined cells into ``4R`` children.

    Child-corner φ is gathered from the parents' corner φ + interior probes via
    ``_CHILD_CORNER_RECIPE`` (no field re-evaluation), with per-child support
    sentinel masking.
    """
    sentinel = _exterior_sentinel(spec)
    x0, x1, y0, y1 = bounds[:, 0], bounds[:, 1], bounds[:, 2], bounds[:, 3]
    xm = 0.5 * (x0 + x1)
    ym = 0.5 * (y0 + y1)
    child_bounds = [
        np.stack([x0, xm, y0, ym], axis=1),  # child 0
        np.stack([xm, x1, y0, ym], axis=1),  # child 1
        np.stack([xm, x1, ym, y1], axis=1),  # child 2
        np.stack([x0, xm, ym, y1], axis=1),  # child 3
    ]
    child_corner_phi: dict[str, list[np.ndarray]] = {lbl: [] for lbl in labels}
    for lbl in labels:
        c = corners[lbl]            # (R, 4)
        p = probes[lbl]             # (R, n_probes)
        box = supports.get(lbl) if supports is not None else None
        for k in range(4):
            cols = [
                (c[:, idx] if kind == "c" else p[:, idx])
                for kind, idx in _CHILD_CORNER_RECIPE[k]
            ]
            vals = np.stack(cols, axis=1)  # (R, 4)
            if box is not None:
                bx0, bx1, by0, by1 = box
                cb = child_bounds[k]
                overlap = ~(
                    (cb[:, 1] <= bx0)
                    | (cb[:, 0] >= bx1)
                    | (cb[:, 3] <= by0)
                    | (cb[:, 2] >= by1)
                )
                vals[~overlap] = sentinel
            child_corner_phi[lbl].append(vals)
    out_bounds = np.concatenate(child_bounds, axis=0)
    out_levels = np.concatenate([levels + 1] * 4, axis=0)
    out_corners = {
        lbl: np.concatenate(child_corner_phi[lbl], axis=0) for lbl in labels
    }
    return out_bounds, out_levels, out_corners
