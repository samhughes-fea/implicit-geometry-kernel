"""The CSG shapes of the kernel: one class per shape, each an ``SDFNode``.

Every shape is a frozen, hashable dataclass that can be evaluated directly
(``evaluate`` on ``(N, 3)`` points, ``z`` ignored by the 2-D shapes), knows its
bounding box (:meth:`Expr.bbox`) and can be shared between trees. Combinators
propagate feature ids (``meta``) the same way as the tagged nodes in
:mod:`implicit_geometry_kernel.sdf.nodes`: a union keeps the id of the nearer child,
an intersection the id of the farther one, and offsets, slabs, shells, rotations and
subtractions keep the id of their base.

There is no separate compile step: a shape is its own evaluator.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

from implicit_geometry_kernel.sdf.field import NO_FEATURE, EvalCache, FieldSample, SDFNode
from implicit_geometry_kernel.sdf.primitives import (
    sdf_box,
    sdf_capsule,
    sdf_circle,
    sdf_half_plane,
    sdf_oriented_box,
    sdf_polygon,
    sdf_rounded_box,
)

BBox = Tuple[float, float, float, float]


@dataclass(frozen=True)
class Expr(SDFNode):
    """Base class of the CSG shapes."""

    def bbox(self) -> Optional[BBox]:
        """Conservative ``(x_min, x_max, y_min, y_max)`` of the region, or ``None``."""
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


def _pick(meta_a, meta_b, take_a) -> np.ndarray:
    return np.where(take_a, meta_a, meta_b).astype(np.int32)


@dataclass(frozen=True)
class Union(Expr):
    """``min``; the feature id of the nearer child."""

    left: Expr
    right: Expr

    def children(self):
        return (self.left, self.right)

    def _evaluate(self, points, cache):
        a, b = self.left.evaluate(points, cache), self.right.evaluate(points, cache)
        return FieldSample(np.minimum(a.phi, b.phi), _pick(a.meta, b.meta, a.phi <= b.phi))


@dataclass(frozen=True)
class Intersect(Expr):
    """``max``; the feature id of the farther child."""

    left: Expr
    right: Expr

    def children(self):
        return (self.left, self.right)

    def _evaluate(self, points, cache):
        a, b = self.left.evaluate(points, cache), self.right.evaluate(points, cache)
        return FieldSample(np.maximum(a.phi, b.phi), _pick(a.meta, b.meta, a.phi >= b.phi))


@dataclass(frozen=True)
class Subtract(Expr):
    """``max(base, -cutter)``; the feature id of the base."""

    base: Expr
    cutter: Expr

    def children(self):
        return (self.base, self.cutter)

    def _evaluate(self, points, cache):
        b, c = self.base.evaluate(points, cache), self.cutter.evaluate(points, cache)
        return FieldSample(np.maximum(b.phi, -c.phi), b.meta)


@dataclass(frozen=True)
class Offset(Expr):
    """``base - amount`` (a positive amount grows the region)."""

    base: Expr
    amount: float

    def children(self):
        return (self.base,)

    def _evaluate(self, points, cache):
        s = self.base.evaluate(points, cache)
        return FieldSample(s.phi - float(self.amount), s.meta)


@dataclass(frozen=True)
class Rotate(Expr):
    """``base`` turned by ``angle`` (radians, anticlockwise) about ``(cx, cy)``."""

    base: Expr
    angle: float
    cx: float = 0.0
    cy: float = 0.0

    def children(self):
        return (self.base,)

    def _evaluate(self, points, cache):
        c, s = np.cos(-float(self.angle)), np.sin(-float(self.angle))
        x, y = points[:, 0] - float(self.cx), points[:, 1] - float(self.cy)
        q = np.column_stack([c * x - s * y + float(self.cx), s * x + c * y + float(self.cy), points[:, 2]])
        return self.base.evaluate(q, EvalCache(q))


@dataclass(frozen=True)
class Shell(Expr):
    """Shell band: ``|base| - thickness / 2``."""

    base: Expr
    thickness: float

    def children(self):
        return (self.base,)

    def _evaluate(self, points, cache):
        s = self.base.evaluate(points, cache)
        return FieldSample(np.abs(s.phi) - float(self.thickness) / 2.0, s.meta)


@dataclass(frozen=True)
class Slab(Expr):
    """Wall slab between base=0 and base=-depth: ``max(base, -base - depth)``."""

    base: Expr
    depth: float

    def children(self):
        return (self.base,)

    def _evaluate(self, points, cache):
        s = self.base.evaluate(points, cache)
        return FieldSample(np.maximum(s.phi, -s.phi - float(self.depth)), s.meta)


def union_expr(exprs: Sequence[Expr]) -> Expr:
    """Balanced n-ary union tree."""
    items = list(exprs)
    if not items:
        raise ValueError("union_expr requires at least one expression.")
    while len(items) > 1:
        nxt: list[Expr] = []
        for i in range(0, len(items), 2):
            if i + 1 < len(items):
                nxt.append(Union(items[i], items[i + 1]))
            else:
                nxt.append(items[i])
        items = nxt
    return items[0]


def intersect_expr(exprs: Sequence[Expr]) -> Expr:
    """Balanced n-ary intersection tree."""
    items = list(exprs)
    if not items:
        raise ValueError("intersect_expr requires at least one expression.")
    while len(items) > 1:
        nxt: list[Expr] = []
        for i in range(0, len(items), 2):
            if i + 1 < len(items):
                nxt.append(Intersect(items[i], items[i + 1]))
            else:
                nxt.append(items[i])
        items = nxt
    return items[0]


def subtract_expr(base: Expr, cutters: Sequence[Expr]) -> Expr:
    """Subtract one or more cutters from a base expression."""
    out = base
    for cutter in cutters:
        out = Subtract(out, cutter)
    return out


def _bbox_union(a: Optional[BBox], b: Optional[BBox]) -> Optional[BBox]:
    if a is None or b is None:
        return None
    return (min(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), max(a[3], b[3]))


def _bbox_intersect(a: Optional[BBox], b: Optional[BBox]) -> Optional[BBox]:
    # Region(A ∩ B) ⊆ Region(A) and ⊆ Region(B), so the box intersection is a
    # valid (and tighter) bound. If only one child is bounded, it still bounds
    # the result.
    if a is None:
        return b
    if b is None:
        return a
    lo_x, hi_x = max(a[0], b[0]), min(a[1], b[1])
    lo_y, hi_y = max(a[2], b[2]), min(a[3], b[3])
    if lo_x > hi_x or lo_y > hi_y:
        # Disjoint boxes ⇒ the intersection region is empty; return a degenerate
        # box so downstream culling drops the (empty) component.
        return (lo_x, lo_x, lo_y, lo_y)
    return (lo_x, hi_x, lo_y, hi_y)


def _bbox_expand(box: Optional[BBox], margin: float) -> Optional[BBox]:
    if box is None:
        return None
    m = float(margin)
    return (box[0] - m, box[1] + m, box[2] - m, box[3] + m)


def _bbox_rotate(box: Optional[BBox], angle: float, cx: float, cy: float) -> Optional[BBox]:
    if box is None:
        return None
    cos_a = math.cos(float(angle))
    sin_a = math.sin(float(angle))
    cxf, cyf = float(cx), float(cy)
    xs: list[float] = []
    ys: list[float] = []
    for px in (box[0], box[1]):
        for py in (box[2], box[3]):
            xr, yr = px - cxf, py - cyf
            xs.append(cxf + cos_a * xr - sin_a * yr)
            ys.append(cyf + sin_a * xr + cos_a * yr)
    return (min(xs), max(xs), min(ys), max(ys))


def expr_bbox(expr: Expr) -> Optional[BBox]:
    """Conservative axis-aligned bbox ``(x_min, x_max, y_min, y_max)`` of the
    region ``expr <= 0`` (the component support), or ``None`` if the region is
    unbounded (``HalfPlane``) or opaque (``CallableField``).

    The box is a guaranteed superset of the true region, so it is safe for
    support-based cell culling: a cell is never dropped if it could contain part
    of the component, and only strictly-exterior cells are culled. Derived
    analytically from the IR graph — no SDF evaluations.
    """
    if isinstance(expr, Circle):
        return (expr.cx - expr.r, expr.cx + expr.r, expr.cy - expr.r, expr.cy + expr.r)
    if isinstance(expr, Box):
        return (
            expr.cx - expr.half_w,
            expr.cx + expr.half_w,
            expr.cy - expr.half_h,
            expr.cy + expr.half_h,
        )
    if isinstance(expr, RoundedBox):
        return (
            expr.cx - expr.half_w - expr.r,
            expr.cx + expr.half_w + expr.r,
            expr.cy - expr.half_h - expr.r,
            expr.cy + expr.half_h + expr.r,
        )
    if isinstance(expr, OrientedBox):
        base = (
            expr.cx - expr.half_w,
            expr.cx + expr.half_w,
            expr.cy - expr.half_h,
            expr.cy + expr.half_h,
        )
        return _bbox_rotate(base, expr.angle, expr.cx, expr.cy)
    if isinstance(expr, Capsule):
        return (
            min(expr.ax, expr.bx) - expr.r,
            max(expr.ax, expr.bx) + expr.r,
            min(expr.ay, expr.by) - expr.r,
            max(expr.ay, expr.by) + expr.r,
        )
    if isinstance(expr, Polygon):
        verts = np.asarray(expr.vertices, dtype=float)
        return (
            float(verts[:, 0].min()),
            float(verts[:, 0].max()),
            float(verts[:, 1].min()),
            float(verts[:, 1].max()),
        )
    if isinstance(expr, (HalfPlane, CallableField)):
        return None
    if isinstance(expr, Union):
        return _bbox_union(expr_bbox(expr.left), expr_bbox(expr.right))
    if isinstance(expr, Intersect):
        return _bbox_intersect(expr_bbox(expr.left), expr_bbox(expr.right))
    if isinstance(expr, Subtract):
        # Region ⊆ Region(base); the cutter only removes interior.
        return expr_bbox(expr.base)
    if isinstance(expr, Offset):
        # {base <= amount}. Positive amount dilates the region (valid only when
        # base is a true distance field, which holds for the primitive leaves
        # this wraps); negative amount erodes ⇒ base box already bounds it.
        return _bbox_expand(expr_bbox(expr.base), max(float(expr.amount), 0.0))
    if isinstance(expr, Shell):
        # |base| - t/2 <= 0 ⇒ band of half-width t/2 around the base boundary.
        return _bbox_expand(expr_bbox(expr.base), abs(float(expr.thickness)) / 2.0)
    if isinstance(expr, Slab):
        # max(base, -base - d) <= 0 ⇒ {-d <= base <= 0} ⊆ Region(base).
        return expr_bbox(expr.base)
    if isinstance(expr, Rotate):
        return _bbox_rotate(expr_bbox(expr.base), expr.angle, expr.cx, expr.cy)
    return None


def eval_expr(expr: Expr, x, y):
    """Evaluate an expression on arbitrary coordinate arrays (no memoization)."""
    return eval_expr_batch(expr, x, y, memo=None)


def eval_expr_batch(
    expr: Expr,
    x,
    y,
    memo: Dict[Expr, np.ndarray] | None = None,
) -> np.ndarray:
    """``phi`` of ``expr`` at coordinate arrays ``x`` and ``y`` (same shape out).

    With ``memo``, every subtree's values are stored by subtree, so a shape shared by
    several trees (one airfoil under every wall) is evaluated once per point set.
    """
    xa, ya = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    if memo is None:
        return _sample(expr, xa, ya)
    return _eval_memo(expr, xa, ya, memo)


def _sample(expr: Expr, xa: np.ndarray, ya: np.ndarray) -> np.ndarray:
    xy = np.column_stack([xa.ravel(), ya.ravel()])
    return np.asarray(expr.evaluate_plane(xy, 0.0).phi, dtype=float).reshape(xa.shape)


def _eval_memo(e: Expr, xs: np.ndarray, ys: np.ndarray, memo: Dict[Expr, np.ndarray]) -> np.ndarray:
    hit = memo.get(e)
    if hit is not None:
        return hit
    if isinstance(e, Union):
        out = np.minimum(_eval_memo(e.left, xs, ys, memo), _eval_memo(e.right, xs, ys, memo))
    elif isinstance(e, Intersect):
        out = np.maximum(_eval_memo(e.left, xs, ys, memo), _eval_memo(e.right, xs, ys, memo))
    elif isinstance(e, Subtract):
        out = np.maximum(_eval_memo(e.base, xs, ys, memo), -_eval_memo(e.cutter, xs, ys, memo))
    elif isinstance(e, Offset):
        out = _eval_memo(e.base, xs, ys, memo) - float(e.amount)
    elif isinstance(e, Shell):
        out = np.abs(_eval_memo(e.base, xs, ys, memo)) - float(e.thickness) / 2.0
    elif isinstance(e, Slab):
        base = _eval_memo(e.base, xs, ys, memo)
        out = np.maximum(base, -base - float(e.depth))
    elif isinstance(e, Rotate):
        c, s = np.cos(-float(e.angle)), np.sin(-float(e.angle))
        xr, yr = xs - float(e.cx), ys - float(e.cy)
        out = _eval_memo(e.base, c * xr - s * yr + float(e.cx), s * xr + c * yr + float(e.cy), memo)
    else:
        out = _sample(e, xs, ys)
    out = np.asarray(out, dtype=float)
    memo[e] = out
    return out


def compile_expr(expr: Expr, grid, cache: Dict[Expr, np.ndarray] | None = None) -> np.ndarray:
    """Evaluate ``expr`` on a grid with subtree memoization."""
    return eval_expr_batch(expr, grid.X, grid.Y, {} if cache is None else cache)


__all__ = [
    "Expr",
    "Circle",
    "Box",
    "RoundedBox",
    "OrientedBox",
    "Capsule",
    "HalfPlane",
    "Polygon",
    "CallableField",
    "Union",
    "Intersect",
    "Subtract",
    "Offset",
    "Rotate",
    "Shell",
    "Slab",
    "union_expr",
    "intersect_expr",
    "subtract_expr",
    "expr_bbox",
    "BBox",
    "eval_expr",
    "eval_expr_batch",
    "compile_expr",
]
