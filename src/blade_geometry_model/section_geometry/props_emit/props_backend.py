"""Property backend configuration for section_geometry."""

from __future__ import annotations

import os
from typing import Literal

import numpy as np

from blade_geometry_model.area_2d.spec import QuadtreeSpec

PropsBackend = Literal["grid", "quadtree", "hybrid"]


def _rotation_invariant_min_extent(verts: np.ndarray) -> float:
    """Minor-axis extent of planar points (rotation-invariant, ~max thickness)."""
    v = np.asarray(verts, dtype=float).reshape(-1, 2)
    if len(v) < 3:
        return float(min(np.ptp(v[:, 0]), np.ptp(v[:, 1])))
    c = v - v.mean(axis=0)
    _, vecs = np.linalg.eigh(np.cov(c.T))
    proj = c @ vecs
    return float(np.ptp(proj, axis=0).min())


def props_backend_from_env(default: PropsBackend = "hybrid") -> PropsBackend:
    raw = os.getenv("SECTION_GEOMETRY_PROPS_BACKEND", default).strip().lower()
    if raw in {"grid", "uniform", "uniform_grid"}:
        return "grid"
    if raw in {"quadtree", "qt", "adaptive"}:
        return "quadtree"
    if raw in {"hybrid", "constructive", "dual"}:
        return "hybrid"
    return default


def quadtree_spec_for_airfoil(
    airfoil_sdf,
    *,
    chord: float,
    padding: float = 0.1,
    max_depth: int | None = None,
    min_cell_size: float | None = None,
    refine_band: float | None = None,
    max_cells: int | None = None,
    root_cell_size: float | None = None,
) -> QuadtreeSpec:
    """Default ``QuadtreeSpec`` sized to an airfoil bounding box."""
    verts = np.asarray(getattr(airfoil_sdf, "vertices", airfoil_sdf), dtype=float)
    depth = int(max_depth if max_depth is not None else os.getenv("SECTION_GEOMETRY_QT_DEPTH", "7"))
    mcs = float(min_cell_size if min_cell_size is not None else 0.0015 * float(chord))
    band = float(refine_band if refine_band is not None else 0.012 * float(chord))
    mc = int(max_cells if max_cells is not None else os.getenv("SECTION_GEOMETRY_QT_MAX_CELLS", "150000"))
    band_scale = float(os.getenv("SECTION_GEOMETRY_QT_BAND_SCALE", "1.5"))
    if root_cell_size is None:
        minor = _rotation_invariant_min_extent(verts)
        root_cell_size = max(minor * (1.0 + 2.0 * padding), 4.0 * mcs)
    return QuadtreeSpec.from_points(
        verts,
        padding=padding,
        max_depth=depth,
        min_cell_size=mcs,
        refine_band=band,
        refine_band_scale=band_scale,
        refine_probes="minimal",
        max_cells=mc,
        square_cells=True,
        root_cell_size=float(root_cell_size),
    )


def quadtree_spec_to_dict(spec: QuadtreeSpec) -> dict[str, float | int | bool | None]:
    return {
        "x_min": float(spec.x_min),
        "x_max": float(spec.x_max),
        "y_min": float(spec.y_min),
        "y_max": float(spec.y_max),
        "max_depth": int(spec.max_depth),
        "min_cell_size": float(spec.min_cell_size),
        "refine_band": float(spec.refine_band),
        "refine_band_scale": float(spec.refine_band_scale),
        "refine_probes": str(spec.refine_probes),
        "max_cells": int(spec.max_cells) if spec.max_cells is not None else None,
        "square_cells": bool(spec.square_cells),
        "root_cell_size": float(spec.root_cell_size) if spec.root_cell_size is not None else None,
    }
