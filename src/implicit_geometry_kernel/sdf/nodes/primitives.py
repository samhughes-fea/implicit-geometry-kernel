"""CSG shapes, part 1: the base class and the 2-D primitive leaves.

Every shape is a frozen, hashable dataclass and an ``SDFNode``: it evaluates itself
(``evaluate`` on ``(N, 3)`` points, ``z`` ignored), knows its bounding box
(:meth:`Expr.bbox`) and can be shared between trees. Combinators are in
:mod:`.csg`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Tuple

import numpy as np

from implicit_geometry_kernel.sdf.kernels.primitives import (
    sdf_box,
    sdf_capsule,
    sdf_circle,
    sdf_half_plane,
    sdf_oriented_box,
    sdf_polygon,
    sdf_rounded_box,
)
from implicit_geometry_kernel.sdf.nodes.base import NO_FEATURE, EvalCache, FieldSample, SDFNode

BBox = Tuple[float, float, float, float]


@dataclass(frozen=True)
class Expr(SDFNode):
    """Base class of the CSG shapes."""

    def bbox(self) -> Optional[BBox]:
        """Conservative ``(x_min, x_max, y_min, y_max)`` of the region, or ``None``."""
        from implicit_geometry_kernel.sdf.nodes.bbox import expr_bbox

        return expr_bbox(self)


@dataclass(frozen=True)
class _Leaf(Expr):
    """A 2-D primitive: ``phi`` from ``(x, y)``, no feature id."""

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        phi = np.asarray(self._phi(points[:, 0], points[:, 1]), dtype=np.float64)
        return FieldSample(phi, np.full(points.shape[0], NO_FEATURE, dtype=np.int32))

    def _phi(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        raise NotImplementedError


@dataclass(frozen=True)
class Circle(_Leaf):
    cx: float = 0.0
    cy: float = 0.0
    r: float = 1.0

    def _phi(self, x, y):
        return sdf_circle(x, y, self.cx, self.cy, self.r)


@dataclass(frozen=True)
class Box(_Leaf):
    cx: float = 0.0
    cy: float = 0.0
    half_w: float = 1.0
    half_h: float = 0.5

    def _phi(self, x, y):
        return sdf_box(x, y, cx=self.cx, cy=self.cy, half_w=self.half_w, half_h=self.half_h)


@dataclass(frozen=True)
class RoundedBox(_Leaf):
    cx: float = 0.0
    cy: float = 0.0
    half_w: float = 1.0
    half_h: float = 0.5
    r: float = 0.1

    def _phi(self, x, y):
        return sdf_rounded_box(x, y, cx=self.cx, cy=self.cy, half_w=self.half_w, half_h=self.half_h, r=self.r)


@dataclass(frozen=True)
class OrientedBox(_Leaf):
    cx: float = 0.0
    cy: float = 0.0
    half_w: float = 1.0
    half_h: float = 0.5
    angle: float = 0.0

    def _phi(self, x, y):
        return sdf_oriented_box(x, y, cx=self.cx, cy=self.cy, half_w=self.half_w, half_h=self.half_h,
                                angle=self.angle)


@dataclass(frozen=True)
class Capsule(_Leaf):
    ax: float
    ay: float
    bx: float
    by: float
    r: float

    def _phi(self, x, y):
        return sdf_capsule(x, y, ax=self.ax, ay=self.ay, bx=self.bx, by=self.by, r=self.r)


@dataclass(frozen=True)
class HalfPlane(_Leaf):
    nx: float = 0.0
    ny: float = 1.0
    d: float = 0.0

    def _phi(self, x, y):
        return sdf_half_plane(x, y, nx=self.nx, ny=self.ny, d=self.d)


@dataclass(frozen=True)
class Polygon(_Leaf):
    vertices: Tuple[Tuple[float, float], ...]
    _hash: int = field(init=False, compare=False, repr=False)
    _array: np.ndarray = field(init=False, compare=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "_hash", hash(self.vertices))
        object.__setattr__(self, "_array", np.asarray(self.vertices, dtype=float))

    def __hash__(self) -> int:
        return self._hash

    def _phi(self, x, y):
        return sdf_polygon(x, y, self._array)


@dataclass(frozen=True, eq=False)
class CallableField(_Leaf):
    """Opaque callable leaf: a token plus the callable itself.

    Two leaves are equal only when they share the token *and* the same callable
    object. Comparing the token alone let a memo keyed on one section's callable
    answer for another section whose token happened to match (tokens built from
    ``id()`` repeat once the first object is collected). Holding ``fn`` in the leaf
    keeps its id from being reused while the leaf, or any memo keyed on it, is alive.
    """

    key: str
    fn: Any = field(repr=False)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CallableField):
            return NotImplemented
        return self.key == other.key and self.fn is other.fn

    def __hash__(self) -> int:
        return hash((self.key, id(self.fn)))

    def _phi(self, x, y):
        return self.fn(x, y)


__all__ = [
    "BBox", "Box", "CallableField", "Capsule", "Circle", "Expr", "HalfPlane", "OrientedBox",
    "Polygon", "RoundedBox",
]
