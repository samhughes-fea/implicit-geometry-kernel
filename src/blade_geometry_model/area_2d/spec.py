"""Quadtree area-mesh specification."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class QuadtreeSpec:
    """Configuration for adaptive quadtree area meshing."""

    x_min: float
    x_max: float
    y_min: float
    y_max: float
    max_depth: int = 6
    min_cell_size: float = 0.005
    refine_band: float = 0.03
    refine_band_scale: float = 1.5
    refine_probes: str = "minimal"
    max_cells: int = 20000
    square_cells: bool = True
    root_cell_size: float | None = None

    @classmethod
    def from_points(
        cls,
        pts: np.ndarray,
        *,
        padding: float = 0.05,
        **kwargs: object,
    ) -> QuadtreeSpec:
        """Bounding box from planar points ``(N, 2)`` in metres with relative padding."""
        arr = np.asarray(pts, dtype=float).reshape(-1, 2)
        xmin, ymin = arr.min(axis=0)
        xmax, ymax = arr.max(axis=0)
        span_x = float(xmax - xmin) + 1e-12
        span_y = float(ymax - ymin) + 1e-12
        pad_x = padding * span_x
        pad_y = padding * span_y
        kw = {k: v for k, v in kwargs.items() if k in cls.__dataclass_fields__}
        return cls(
            x_min=float(xmin - pad_x),
            x_max=float(xmax + pad_x),
            y_min=float(ymin - pad_y),
            y_max=float(ymax + pad_y),
            **kw,  # type: ignore[arg-type]
        )


# Backward-compatible alias (root_blend_model naming).
QuadTreeSpec = QuadtreeSpec
