"""CSG shapes, part 2: combinators and evaluation helpers.

Combinators propagate feature ids (``meta``): a union keeps the id of the nearer
child, an intersection the id of the farther one, and offsets, slabs, shells,
rotations and subtractions keep the id of their base.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Sequence

import numpy as np

from implicit_geometry_kernel.sdf.nodes.base import EvalCache, FieldSample
from implicit_geometry_kernel.sdf.nodes.primitives import Expr


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
    "Intersect", "Offset", "Rotate", "Shell", "Slab", "Subtract", "Union", "compile_expr",
    "eval_expr", "eval_expr_batch", "intersect_expr", "subtract_expr", "union_expr",
]
