"""Adaptive quadtree mesh for scalar-field zero-contour extraction."""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field, replace

from typing import Any, Callable

import numpy as np

from .contour import ensure_closed, smooth_contour
from .spec import QuadtreeSpec

SupportBox = tuple[float, float, float, float]

# Shared read-only "missing label" sentinel for _phi_by_label_from_leaves. Hoisted
# so it is not re-allocated as a dict.get default on every leaf x label (np.array
# copies it into the result, so sharing the object is safe and bit-identical).
_NAN4 = np.full(4, np.nan)


@dataclass
class QuadCell:
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    level: int = 0
    phi_corners: np.ndarray | None = None
    phi_fields: dict[str, np.ndarray] | None = None
    children: list[QuadCell] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0

    @property
    def width(self) -> float:
        return self.x_max - self.x_min

    @property
    def height(self) -> float:
        return self.y_max - self.y_min

    def corners(self) -> np.ndarray:
        return np.array(
            [
                [self.x_min, self.y_min],
                [self.x_max, self.y_min],
                [self.x_max, self.y_max],
                [self.x_min, self.y_max],
            ],
            dtype=float,
        )


@dataclass
class QuadTreeMesh:
    spec: QuadtreeSpec
    roots: list[QuadCell]
    refine_stats: dict[str, int] = field(
        default_factory=lambda: {
            "sign_change": 0,
            "near_interface": 0,
            "subdivisions": 0,
        }
    )

    def leaf_cells(self) -> list[QuadCell]:
        cached = getattr(self, "_cached_leaves", None)
        if cached is not None:
            return cached
        if getattr(self, "_is_soa", False):
            return _materialize_soa_leaves(self)
        return _refresh_leaf_cache(self)


def _materialize_soa_leaves(mesh: QuadTreeMesh) -> list[QuadCell]:
    """Lazily build ``QuadCell`` objects from a struct-of-arrays mesh.

    SoA meshes (from ``build_shared_quadtree``) store leaves as ``_leaf_bounds`` /
    ``_leaf_levels`` / ``phi_by_label`` arrays; the hot property/plot paths read
    those directly. Object-iterating consumers (contour stitching, the rare
    hollow-cell sampler, examples) call ``leaf_cells()``, which materialises here
    on demand and caches the result.
    """
    bounds = mesh._leaf_bounds  # type: ignore[attr-defined]
    levels = mesh._leaf_levels  # type: ignore[attr-defined]
    phi = getattr(mesh, "phi_by_label", {}) or {}
    labels = list(phi.keys())
    single = len(labels) == 1
    leaves: list[QuadCell] = []
    for i in range(bounds.shape[0]):
        cell = QuadCell(
            float(bounds[i, 0]),
            float(bounds[i, 1]),
            float(bounds[i, 2]),
            float(bounds[i, 3]),
            level=int(levels[i]),
        )
        cell.phi_fields = {lbl: phi[lbl][i] for lbl in labels}
        if single:
            cell.phi_corners = next(iter(cell.phi_fields.values()))
        leaves.append(cell)
    mesh._cached_leaves = leaves  # type: ignore[attr-defined]
    return leaves


def _invalidate_leaf_cache(mesh: QuadTreeMesh) -> None:
    mesh._cached_leaves = None  # type: ignore[attr-defined]
    mesh._leaf_bounds = None  # type: ignore[attr-defined]


def _refresh_leaf_cache(mesh: QuadTreeMesh) -> list[QuadCell]:
    leaves: list[QuadCell] = []
    stack = list(mesh.roots)
    while stack:
        cell = stack.pop()
        if cell.is_leaf:
            leaves.append(cell)
        else:
            stack.extend(cell.children)
    mesh._cached_leaves = leaves  # type: ignore[attr-defined]
    if leaves:
        mesh._leaf_bounds = np.array(  # type: ignore[attr-defined]
            [[c.x_min, c.x_max, c.y_min, c.y_max] for c in leaves],
            dtype=float,
        )
    else:
        mesh._leaf_bounds = np.empty((0, 4), dtype=float)  # type: ignore[attr-defined]
    return leaves


def leaf_bounds(mesh: QuadTreeMesh) -> np.ndarray:
    """Cached leaf bounds ``(n_leaves, 4)`` as x_min, x_max, y_min, y_max."""
    bounds = getattr(mesh, "_leaf_bounds", None)
    if bounds is not None:
        return bounds
    _refresh_leaf_cache(mesh)
    return getattr(mesh, "_leaf_bounds", np.empty((0, 4), dtype=float))


def _cells_bounds_array(cells: list[QuadCell]) -> np.ndarray:
    if not cells:
        return np.empty((0, 4), dtype=float)
    return np.array(
        [[c.x_min, c.x_max, c.y_min, c.y_max] for c in cells],
        dtype=float,
    )


def _sample_cell_phi(cell: QuadCell, phi_sampler) -> np.ndarray:
    return np.asarray(phi_sampler(cell.corners()), dtype=float)


def _domain_span(spec: QuadtreeSpec) -> float:
    return max(float(spec.x_max - spec.x_min), float(spec.y_max - spec.y_min), 1e-12)


def _exterior_sentinel(spec: QuadtreeSpec) -> float:
    return 2.0 * _domain_span(spec)


def _cell_diagonal(cell: QuadCell) -> float:
    return math.hypot(cell.width, cell.height)


def _cell_overlaps_support(cell: QuadCell, support: SupportBox) -> bool:
    xmin, xmax, ymin, ymax = support
    return not (
        cell.x_max <= xmin or cell.x_min >= xmax or cell.y_max <= ymin or cell.y_min >= ymax
    )


def _n_interior_probes(spec: QuadtreeSpec) -> int:
    if str(spec.refine_probes).lower() == "full":
        return 13
    return 5


def _cell_interior_points(cell: QuadCell, spec: QuadtreeSpec) -> np.ndarray:
    xm = 0.5 * (cell.x_min + cell.x_max)
    ym = 0.5 * (cell.y_min + cell.y_max)
    if str(spec.refine_probes).lower() == "full":
        xq1 = cell.x_min + 0.25 * cell.width
        xq3 = cell.x_min + 0.75 * cell.width
        yq1 = cell.y_min + 0.25 * cell.height
        yq3 = cell.y_min + 0.75 * cell.height
        return np.array(
            [
                [xm, ym],
                [xm, cell.y_min],
                [cell.x_max, ym],
                [xm, cell.y_max],
                [cell.x_min, ym],
                [xq1, cell.y_min],
                [xq3, cell.y_min],
                [cell.x_max, yq1],
                [cell.x_max, yq3],
                [xq1, cell.y_max],
                [xq3, cell.y_max],
                [cell.x_min, yq1],
                [cell.x_min, yq3],
            ],
            dtype=float,
        )
    return np.array(
        [
            [xm, ym],
            [xm, cell.y_min],
            [cell.x_max, ym],
            [xm, cell.y_max],
            [cell.x_min, ym],
        ],
        dtype=float,
    )


def _subdivide_cell(cell: QuadCell, phi_sampler) -> list[QuadCell]:
    xm = 0.5 * (cell.x_min + cell.x_max)
    ym = 0.5 * (cell.y_min + cell.y_max)
    children = [
        QuadCell(cell.x_min, xm, cell.y_min, ym, level=cell.level + 1),
        QuadCell(xm, cell.x_max, cell.y_min, ym, level=cell.level + 1),
        QuadCell(xm, cell.x_max, ym, cell.y_max, level=cell.level + 1),
        QuadCell(cell.x_min, xm, ym, cell.y_max, level=cell.level + 1),
    ]
    all_corners = np.vstack([c.corners() for c in children])
    all_phi = np.asarray(phi_sampler(all_corners), dtype=float)
    for i, child in enumerate(children):
        child.phi_corners = all_phi[i * 4 : (i + 1) * 4]
    cell.children = children
    return children


def build_root_quadtree(spec: QuadtreeSpec, phi_sampler) -> QuadTreeMesh:
    """Create a single-cell root mesh and sample corner phi values."""
    root = QuadCell(spec.x_min, spec.x_max, spec.y_min, spec.y_max, level=0)
    root.phi_corners = _sample_cell_phi(root, phi_sampler)
    return QuadTreeMesh(spec=spec, roots=[root])


def _square_forest_domain(spec: QuadtreeSpec) -> tuple[float, float, float, float, float, int, int]:
    """Return expanded square-tiled domain and grid counts for a root forest."""
    width = float(spec.x_max - spec.x_min)
    height = float(spec.y_max - spec.y_min)
    s = float(spec.root_cell_size) if spec.root_cell_size is not None else min(width, height)
    if s <= 0.0:
        raise ValueError("root_cell_size must be positive.")
    n_x = max(1, int(math.ceil(width / s)))
    n_y = max(1, int(math.ceil(height / s)))
    cx = 0.5 * (spec.x_min + spec.x_max)
    cy = 0.5 * (spec.y_min + spec.y_max)
    total_w = n_x * s
    total_h = n_y * s
    x_min = cx - 0.5 * total_w
    x_max = cx + 0.5 * total_w
    y_min = cy - 0.5 * total_h
    y_max = cy + 0.5 * total_h
    return x_min, x_max, y_min, y_max, s, n_x, n_y


def build_square_root_forest(spec: QuadtreeSpec, samplers: dict[str, Any]) -> QuadTreeMesh:
    """Tile the domain with equal square root cells (forest of square quadtrees)."""
    if not samplers:
        raise ValueError("samplers must be non-empty.")
    x_min, x_max, y_min, y_max, s, n_x, n_y = _square_forest_domain(spec)
    forest_spec = replace(
        spec,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        root_cell_size=s,
    )
    roots: list[QuadCell] = []
    for j in range(n_y):
        for i in range(n_x):
            roots.append(
                QuadCell(
                    x_min + i * s,
                    x_min + (i + 1) * s,
                    y_min + j * s,
                    y_min + (j + 1) * s,
                    level=0,
                )
            )
    all_corners = np.vstack([c.corners() for c in roots])
    phi_by_root_label: dict[str, np.ndarray] = {}
    for lbl, sampler in samplers.items():
        phi_by_root_label[str(lbl)] = np.asarray(sampler(all_corners), dtype=float)
    for k, cell in enumerate(roots):
        cell.phi_fields = {
            lbl: vals[k * 4 : (k + 1) * 4] for lbl, vals in phi_by_root_label.items()
        }
        if len(samplers) == 1:
            cell.phi_corners = next(iter(cell.phi_fields.values()))
    return QuadTreeMesh(spec=forest_spec, roots=roots)


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
                    for i in active_rows:
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
        for i in notneed_idx:
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
        vals = np.empty(len(leaves), dtype=float)
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


BatchEval = Callable[[np.ndarray], dict[str, np.ndarray]]


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
    sentinel masking — the array-native form of ``_subdivide_refined_reuse``.
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


def _refine_soa(
    spec: QuadtreeSpec,
    labels: list[str],
    batch_eval: BatchEval,
    supports: dict[str, SupportBox] | None,
    frontier: tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]],
) -> QuadTreeMesh:
    """Struct-of-arrays adaptive refinement (no per-cell ``QuadCell`` objects).

    Mirrors the ``refine_quadtree_multi`` decision (corner cache + interior probe
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
                all_phi, band_ceiling, band_scale, cell_diags
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


def _interpolate_zero_crossing(
    p1: np.ndarray,
    p2: np.ndarray,
    v1: float,
    v2: float,
) -> np.ndarray:
    if abs(v2 - v1) < 1e-15:
        return 0.5 * (p1 + p2)
    t = -v1 / (v2 - v1)
    return p1 + t * (p2 - p1)


def _connect_segments_to_contour(
    segments: list[tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    if not segments:
        raise ValueError("No zero contour found in quadtree leaf cells.")
    chains: list[list[np.ndarray]] = []
    tol = 1e-9
    for a, b in segments:
        placed = False
        for chain in chains:
            if np.allclose(chain[-1], a, atol=tol, rtol=0.0):
                chain.append(b)
                placed = True
                break
            if np.allclose(chain[-1], b, atol=tol, rtol=0.0):
                chain.append(a)
                placed = True
                break
            if np.allclose(chain[0], b, atol=tol, rtol=0.0):
                chain.insert(0, a)
                placed = True
                break
            if np.allclose(chain[0], a, atol=tol, rtol=0.0):
                chain.insert(0, b)
                placed = True
                break
        if not placed:
            chains.append([a, b])
    merged = True
    while merged and len(chains) > 1:
        merged = False
        for i in range(len(chains)):
            if merged:
                break
            for j in range(i + 1, len(chains)):
                c1, c2 = chains[i], chains[j]
                if np.allclose(c1[-1], c2[0], atol=tol, rtol=0.0):
                    chains[i] = c1 + c2[1:]
                elif np.allclose(c1[-1], c2[-1], atol=tol, rtol=0.0):
                    chains[i] = c1 + c2[-2::-1]
                elif np.allclose(c1[0], c2[-1], atol=tol, rtol=0.0):
                    chains[i] = c2[:-1] + c1
                elif np.allclose(c1[0], c2[0], atol=tol, rtol=0.0):
                    chains[i] = c2[::-1][:-1] + c1
                else:
                    continue
                del chains[j]
                merged = True
                break
    longest = max(chains, key=len)
    return ensure_closed(np.asarray(longest, dtype=float))


def contour_from_quadtree_marching_squares(
    leaves: list[QuadCell],
    passes: int = 2,
) -> np.ndarray:
    """Extract the zero-level contour from quadtree leaf cells."""
    edge_pairs = [(0, 1), (1, 2), (2, 3), (3, 0)]
    leaves_with_phi = [c for c in leaves if c.phi_corners is not None]
    if not leaves_with_phi:
        raise ValueError("No leaf cells with phi_corners found.")
    phi_array = np.array([c.phi_corners for c in leaves_with_phi], dtype=float)
    sign_change_mask = (phi_array.min(axis=1) <= 0.0) & (phi_array.max(axis=1) >= 0.0)
    active = [c for c, active in zip(leaves_with_phi, sign_change_mask) if active]
    segments = []
    for cell in active:
        phi = cell.phi_corners
        corners = cell.corners()
        crossings = []
        for i, j in edge_pairs:
            if (phi[i] <= 0.0) != (phi[j] <= 0.0):
                crossings.append(
                    _interpolate_zero_crossing(corners[i], corners[j], phi[i], phi[j])
                )
        if len(crossings) == 2:
            segments.append((crossings[0], crossings[1]))
        elif len(crossings) == 4:
            segments.append((crossings[0], crossings[1]))
            segments.append((crossings[2], crossings[3]))
    contour = _connect_segments_to_contour(segments)
    return smooth_contour(contour, passes=passes)


def summarize_quadtree(mesh: QuadTreeMesh) -> dict[str, object]:
    """Summary statistics for a refined quadtree mesh."""
    if getattr(mesh, "_is_soa", False):
        bounds = mesh._leaf_bounds  # type: ignore[attr-defined]
        levels_arr = mesh._leaf_levels  # type: ignore[attr-defined]
        widths = bounds[:, 1] - bounds[:, 0] if bounds.size else np.empty(0)
        uniq, counts = np.unique(levels_arr, return_counts=True)
        level_counts = {int(k): int(v) for k, v in zip(uniq, counts)}
        return {
            "leaf_count": int(bounds.shape[0]),
            "level_counts": dict(sorted(level_counts.items())),
            "max_level": int(levels_arr.max()) if levels_arr.size else 0,
            "min_cell_width": float(widths.min()) if widths.size else 0.0,
            "max_cell_width": float(widths.max()) if widths.size else 0.0,
            "refine_stats": dict(mesh.refine_stats),
        }
    leaves = mesh.leaf_cells()
    levels: dict[int, int] = {}
    for cell in leaves:
        levels[cell.level] = levels.get(cell.level, 0) + 1
    return {
        "leaf_count": len(leaves),
        "level_counts": dict(sorted(levels.items())),
        "max_level": max(levels) if levels else 0,
        "min_cell_width": min((c.width for c in leaves), default=0.0),
        "max_cell_width": max((c.width for c in leaves), default=0.0),
        "refine_stats": dict(mesh.refine_stats),
    }
