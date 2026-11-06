"""Conservative bounding boxes of CSG shapes (``expr_bbox``)."""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from implicit_geometry_kernel.sdf.nodes.csg import Intersect, Offset, Rotate, Shell, Slab, Subtract, Union
from implicit_geometry_kernel.sdf.nodes.primitives import (
    BBox,
    Box,
    CallableField,
    Capsule,
    Circle,
    Expr,
    HalfPlane,
    OrientedBox,
    Polygon,
    RoundedBox,
)


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


__all__ = ["expr_bbox"]
