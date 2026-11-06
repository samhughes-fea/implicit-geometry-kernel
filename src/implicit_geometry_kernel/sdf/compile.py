"""Compile a frozen CSG ``Expr`` into an ``SDFNode`` tree."""

from __future__ import annotations

from functools import lru_cache
from typing import Dict

import numpy as np

from implicit_geometry_kernel.sdf.ir import (
    Box,
    CallableField,
    Capsule,
    Circle,
    Expr,
    HalfPlane,
    Intersect,
    Offset,
    OrientedBox,
    Polygon,
    Rotate,
    RoundedBox,
    Shell,
    Slab,
    Subtract,
    Union,
)
from implicit_geometry_kernel.sdf.nodes import (
    CallableSDF,
    CircleSDF,
    IntersectNode,
    OffsetNode,
    PolygonSDF,
    RotateNode,
    ShellNode,
    SlabNode,
    SubtractNode,
    UnionNode,
)
from implicit_geometry_kernel.sdf.primitives import (
    sdf_box,
    sdf_capsule,
    sdf_half_plane,
    sdf_oriented_box,
    sdf_rounded_box,
)


@lru_cache(maxsize=2048)
def _compile_cached(expr: Expr):
    return _compile_impl(expr)


def compile_expr_to_node(expr: Expr):
    """One ``SDFNode`` per ``Expr``; a new primitive is one class plus this mapping."""
    try:
        return _compile_cached(expr)
    except TypeError:
        return _compile_impl(expr)


def _compile_impl(expr: Expr):
    if isinstance(expr, Circle):
        return CircleSDF(cx=expr.cx, cy=expr.cy, r=expr.r)
    if isinstance(expr, Polygon):
        return PolygonSDF(np.asarray(expr.vertices, dtype=float))
    if isinstance(expr, Box):
        return CallableSDF(
            lambda x, y, e=expr: sdf_box(x, y, cx=e.cx, cy=e.cy, half_w=e.half_w, half_h=e.half_h)
        )
    if isinstance(expr, RoundedBox):
        return CallableSDF(
            lambda x, y, e=expr: sdf_rounded_box(
                x, y, cx=e.cx, cy=e.cy, half_w=e.half_w, half_h=e.half_h, r=e.r
            )
        )
    if isinstance(expr, OrientedBox):
        return CallableSDF(
            lambda x, y, e=expr: sdf_oriented_box(
                x, y, cx=e.cx, cy=e.cy, half_w=e.half_w, half_h=e.half_h, angle=e.angle
            )
        )
    if isinstance(expr, Capsule):
        return CallableSDF(
            lambda x, y, e=expr: sdf_capsule(x, y, ax=e.ax, ay=e.ay, bx=e.bx, by=e.by, r=e.r)
        )
    if isinstance(expr, HalfPlane):
        return CallableSDF(
            lambda x, y, e=expr: sdf_half_plane(x, y, nx=e.nx, ny=e.ny, d=e.d)
        )
    if isinstance(expr, CallableField):
        return CallableSDF(expr.fn)
    if isinstance(expr, Union):
        return UnionNode(compile_expr_to_node(expr.left), compile_expr_to_node(expr.right))
    if isinstance(expr, Intersect):
        return IntersectNode(compile_expr_to_node(expr.left), compile_expr_to_node(expr.right))
    if isinstance(expr, Subtract):
        return SubtractNode(compile_expr_to_node(expr.base), compile_expr_to_node(expr.cutter))
    if isinstance(expr, Offset):
        return OffsetNode(compile_expr_to_node(expr.base), expr.amount)
    if isinstance(expr, Rotate):
        return RotateNode(compile_expr_to_node(expr.base), expr.angle, expr.cx, expr.cy)
    if isinstance(expr, Shell):
        return ShellNode(compile_expr_to_node(expr.base), expr.thickness)
    if isinstance(expr, Slab):
        return SlabNode(compile_expr_to_node(expr.base), expr.depth)
    raise TypeError(f"Unsupported Expr node: {type(expr).__name__}")


def _sample_node(expr: Expr, xa: np.ndarray, ya: np.ndarray) -> np.ndarray:
    xy = np.column_stack([xa.ravel(), ya.ravel()])
    phi = compile_expr_to_node(expr).evaluate_plane(xy, 0.0).phi
    return np.asarray(phi, dtype=float).reshape(xa.shape)


def _rotate_coords(x, y, angle: float, cx: float, cy: float):
    cos_a = np.cos(-float(angle))
    sin_a = np.sin(-float(angle))
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    xr = xa - float(cx)
    yr = ya - float(cy)
    xp = cos_a * xr - sin_a * yr + float(cx)
    yp = sin_a * xr + cos_a * yr + float(cy)
    return xp, yp


def eval_via_nodes(
    expr: Expr,
    x,
    y,
    memo: Dict[Expr, np.ndarray] | None = None,
) -> np.ndarray:
    """Evaluate ``expr`` via ``compile_expr_to_node``.

    When ``memo`` is supplied, child ``Expr``s are evaluated first and stored
    so shared subtrees (e.g. one polygon in two combinators) run once.
    """
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    xa, ya = np.broadcast_arrays(xa, ya)
    if memo is None:
        return _sample_node(expr, xa, ya)

    def rec(e: Expr, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
        hit = memo.get(e)
        if hit is not None:
            return hit
        if isinstance(e, Union):
            out = np.minimum(rec(e.left, xs, ys), rec(e.right, xs, ys))
        elif isinstance(e, Intersect):
            out = np.maximum(rec(e.left, xs, ys), rec(e.right, xs, ys))
        elif isinstance(e, Subtract):
            out = np.maximum(rec(e.base, xs, ys), -rec(e.cutter, xs, ys))
        elif isinstance(e, Offset):
            out = rec(e.base, xs, ys) - float(e.amount)
        elif isinstance(e, Shell):
            out = np.abs(rec(e.base, xs, ys)) - float(e.thickness) / 2.0
        elif isinstance(e, Slab):
            base = rec(e.base, xs, ys)
            out = np.maximum(base, -base - float(e.depth))
        elif isinstance(e, Rotate):
            xp, yp = _rotate_coords(xs, ys, e.angle, e.cx, e.cy)
            out = rec(e.base, xp, yp)
        else:
            out = _sample_node(e, xs, ys)
        out = np.asarray(out, dtype=float)
        memo[e] = out
        return out

    return rec(expr, xa, ya)


__all__ = ["compile_expr_to_node", "eval_via_nodes"]
