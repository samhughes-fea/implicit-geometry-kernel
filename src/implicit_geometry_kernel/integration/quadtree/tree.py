"""Quadtree cells and meshes: data types, leaf caches and root construction."""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

from typing import Any

import numpy as np

from implicit_geometry_kernel.integration.spec import QuadtreeSpec


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


