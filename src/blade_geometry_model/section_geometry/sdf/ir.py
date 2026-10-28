"""Compiled CSG expression graph for grid- and batch-based SDF evaluation."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

BBox = Tuple[float, float, float, float]

from .primitives import (
    sdf_box,
    sdf_capsule,
    sdf_circle,
    sdf_half_plane,
    sdf_oriented_box,
    sdf_polygon,
    sdf_rounded_box,
)


@dataclass(frozen=True)
class Expr:
    """Base expression type."""


@dataclass(frozen=True)
class Circle(Expr):
    cx: float = 0.0
    cy: float = 0.0
    r: float = 1.0


@dataclass(frozen=True)
class Box(Expr):
    cx: float = 0.0
    cy: float = 0.0
    half_w: float = 1.0
    half_h: float = 0.5


@dataclass(frozen=True)
class RoundedBox(Expr):
    cx: float = 0.0
    cy: float = 0.0
    half_w: float = 1.0
    half_h: float = 0.5
    r: float = 0.1


@dataclass(frozen=True)
class OrientedBox(Expr):
    cx: float = 0.0
    cy: float = 0.0
    half_w: float = 1.0
    half_h: float = 0.5
    angle: float = 0.0


@dataclass(frozen=True)
class Capsule(Expr):
    ax: float
    ay: float
    bx: float
    by: float
    r: float


@dataclass(frozen=True)
class HalfPlane(Expr):
    nx: float = 0.0
    ny: float = 1.0
    d: float = 0.0


@dataclass(frozen=True)
class Polygon(Expr):
    vertices: Tuple[Tuple[float, float], ...]
    _hash: int = field(init=False, compare=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "_hash", hash(self.vertices))

    def __hash__(self) -> int:
        return self._hash


@dataclass(frozen=True)
class CallableField(Expr):
    """Opaque callable leaf keyed by a stable token."""

    key: str
    fn: Any = field(compare=False, hash=False, repr=False)


@dataclass(frozen=True)
class Union(Expr):
    left: Expr
    right: Expr


@dataclass(frozen=True)
class Intersect(Expr):
    left: Expr
    right: Expr


@dataclass(frozen=True)
class Subtract(Expr):
    base: Expr
    cutter: Expr


@dataclass(frozen=True)
class Offset(Expr):
    base: Expr
    amount: float


@dataclass(frozen=True)
class Rotate(Expr):
    base: Expr
    angle: float
    cx: float = 0.0
    cy: float = 0.0


@dataclass(frozen=True)
class Shell(Expr):
    """Shell band: |base| - thickness/2."""

    base: Expr
    thickness: float


@dataclass(frozen=True)
class Slab(Expr):
    """Wall slab between base=0 and base=-depth: max(base, -base - depth)."""

    base: Expr
    depth: float


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


def _rotate_coords(x, y, angle: float, cx: float = 0.0, cy: float = 0.0):
    cos_a = np.cos(-float(angle))
    sin_a = np.sin(-float(angle))
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    xr = xa - float(cx)
    yr = ya - float(cy)
    xp = cos_a * xr - sin_a * yr + float(cx)
    yp = sin_a * xr + cos_a * yr + float(cy)
    return xp, yp


def eval_expr(expr: Expr, x, y):
    """Evaluate an expression on arbitrary coordinate arrays (no memoization)."""
    return eval_expr_batch(expr, x, y, memo=None)


def eval_expr_batch(
    expr: Expr,
    x,
    y,
    memo: Dict[Expr, np.ndarray] | None = None,
) -> np.ndarray:
    """Evaluate an expression with optional subtree memoization keyed by ``Expr``."""
    local = memo is None
    store: Dict[Expr, np.ndarray] = {} if memo is None else memo
    hit = store.get(expr)
    if hit is not None:
        return hit

    if isinstance(expr, Circle):
        out = sdf_circle(x, y, cx=expr.cx, cy=expr.cy, r=expr.r)
    elif isinstance(expr, Box):
        out = sdf_box(x, y, cx=expr.cx, cy=expr.cy, half_w=expr.half_w, half_h=expr.half_h)
    elif isinstance(expr, RoundedBox):
        out = sdf_rounded_box(
            x, y, cx=expr.cx, cy=expr.cy, half_w=expr.half_w, half_h=expr.half_h, r=expr.r
        )
    elif isinstance(expr, OrientedBox):
        out = sdf_oriented_box(
            x,
            y,
            cx=expr.cx,
            cy=expr.cy,
            half_w=expr.half_w,
            half_h=expr.half_h,
            angle=expr.angle,
        )
    elif isinstance(expr, Capsule):
        out = sdf_capsule(x, y, ax=expr.ax, ay=expr.ay, bx=expr.bx, by=expr.by, r=expr.r)
    elif isinstance(expr, HalfPlane):
        out = sdf_half_plane(x, y, nx=expr.nx, ny=expr.ny, d=expr.d)
    elif isinstance(expr, Polygon):
        verts = np.asarray(expr.vertices, dtype=float)
        out = sdf_polygon(x, y, verts)
    elif isinstance(expr, CallableField):
        out = np.asarray(expr.fn(x, y), dtype=float)
    elif isinstance(expr, Union):
        out = np.minimum(
            eval_expr_batch(expr.left, x, y, store),
            eval_expr_batch(expr.right, x, y, store),
        )
    elif isinstance(expr, Intersect):
        out = np.maximum(
            eval_expr_batch(expr.left, x, y, store),
            eval_expr_batch(expr.right, x, y, store),
        )
    elif isinstance(expr, Subtract):
        out = np.maximum(
            eval_expr_batch(expr.base, x, y, store),
            -eval_expr_batch(expr.cutter, x, y, store),
        )
    elif isinstance(expr, Offset):
        out = eval_expr_batch(expr.base, x, y, store) - float(expr.amount)
    elif isinstance(expr, Rotate):
        xp, yp = _rotate_coords(x, y, expr.angle, expr.cx, expr.cy)
        out = eval_expr_batch(expr.base, xp, yp, store)
    elif isinstance(expr, Shell):
        base = eval_expr_batch(expr.base, x, y, store)
        out = np.abs(base) - float(expr.thickness) / 2.0
    elif isinstance(expr, Slab):
        base = eval_expr_batch(expr.base, x, y, store)
        depth = float(expr.depth)
        out = np.maximum(base, -base - depth)
    else:
        raise TypeError(f"Unsupported Expr node: {type(expr).__name__}")

    out = np.asarray(out, dtype=float)
    store[expr] = out
    if local:
        return out
    return out


def compile_expr(expr: Expr, grid, cache: Dict[Expr, np.ndarray] | None = None) -> np.ndarray:
    """Compile/evaluate expression on a grid with subtree memoization."""
    memo = {} if cache is None else cache
    hit = memo.get(expr)
    if hit is not None:
        return hit

    if isinstance(expr, Union):
        out = np.minimum(
            compile_expr(expr.left, grid, memo), compile_expr(expr.right, grid, memo)
        )
    elif isinstance(expr, Intersect):
        out = np.maximum(
            compile_expr(expr.left, grid, memo), compile_expr(expr.right, grid, memo)
        )
    elif isinstance(expr, Subtract):
        out = np.maximum(
            compile_expr(expr.base, grid, memo), -compile_expr(expr.cutter, grid, memo)
        )
    elif isinstance(expr, Offset):
        out = compile_expr(expr.base, grid, memo) - float(expr.amount)
    elif isinstance(expr, Shell):
        base = compile_expr(expr.base, grid, memo)
        out = np.abs(base) - float(expr.thickness) / 2.0
    elif isinstance(expr, Slab):
        base = compile_expr(expr.base, grid, memo)
        depth = float(expr.depth)
        out = np.maximum(base, -base - depth)
    else:
        out = np.asarray(eval_expr_batch(expr, grid.X, grid.Y, memo), dtype=float)
    memo[expr] = out
    return out


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
