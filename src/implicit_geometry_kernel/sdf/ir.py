"""Moved to :mod:`implicit_geometry_kernel.sdf.nodes` (primitives, csg, bbox); alias until 0.6.0."""

from implicit_geometry_kernel.sdf.nodes.bbox import expr_bbox
from implicit_geometry_kernel.sdf.nodes.csg import (
    Intersect,
    Offset,
    Rotate,
    Shell,
    Slab,
    Subtract,
    Union,
    compile_expr,
    eval_expr,
    eval_expr_batch,
    intersect_expr,
    subtract_expr,
    union_expr,
)
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

__all__ = [
    "Expr", "Circle", "Box", "RoundedBox", "OrientedBox", "Capsule", "HalfPlane", "Polygon",
    "CallableField", "Union", "Intersect", "Subtract", "Offset", "Rotate", "Shell", "Slab",
    "union_expr", "intersect_expr", "subtract_expr", "expr_bbox", "BBox", "eval_expr",
    "eval_expr_batch", "compile_expr",
]
