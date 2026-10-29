"""2-D area mesh result types and sampling protocols."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Protocol

import numpy as np

from .quadtree import QuadTreeMesh


class PhiSampler(Protocol):
    """Sample a scalar field at planar points ``pts`` with shape ``(N, 2)``."""

    def __call__(self, pts: np.ndarray) -> np.ndarray: ...


@dataclass
class FieldMesh2D:
    """Quadtree mesh and optional cached zero contour for one scalar field."""

    mesh: QuadTreeMesh
    contour: np.ndarray | None = None
    stats: dict[str, object] = field(default_factory=dict)


@dataclass
class AreaMesh2D:
    """Named 2-D area meshes (one per field label)."""

    fields: Mapping[str, FieldMesh2D]
