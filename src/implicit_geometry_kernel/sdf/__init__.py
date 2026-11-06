"""
implicit_geometry_kernel.sdf
===================
Generic 2-D signed distance field (SDF) and constructive solid geometry (CSG)
toolbox.  No blade domain knowledge — reusable for any implicit geometry task.

Public API
----------
Primitives:
    sdf_circle, sdf_box, sdf_half_plane, sdf_segment, sdf_capsule,
    sdf_rounded_box, sdf_ellipse, sdf_oriented_box, sdf_polygon

CSG operations:
    union, intersect, subtract, offset, shell, blend,
    smooth_union, smooth_intersect, smooth_subtract,
    union_all, intersect_all

Transforms:
    rotate_field, translate_field, scale_field, mirror_field_x, mirror_field_y,
    rotate_points, flapwise_aligned_web_angle, SDFFrame

Grid:
    SDFGrid

Airfoil:
    AirfoilSDF (adapter; prefer Airfoil.naca4(...).to_sdf())

CSG shapes (each an SDFNode; evaluate directly or with eval_expr_batch):
    Expr, Circle, Box, RoundedBox, OrientedBox, Capsule, HalfPlane, Polygon,
    CallableField, Union, Intersect, Subtract, Offset, Rotate, Shell, Slab,
    eval_expr, eval_expr_batch, compile_expr

Feature field (2.5-D):
    FieldSample, FeatureTable, FeatureSpec, NO_FEATURE, SDFNode,
    TagNode, PriorityUnionNode, ComplementNode, MaxNode, LawOffsetNode,
    SpanIntervalNode, LoftedPolygonSDF, GeometrySolver

Protocol:
    SDFCallable
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable

@runtime_checkable
class SDFCallable(Protocol):
    """Callable signed-distance field interface: (x, y) -> ndarray."""

    def __call__(self, x: Any, y: Any) -> Any: ...

from implicit_geometry_kernel.sdf.kernels.primitives import sdf_box, sdf_capsule, sdf_circle, sdf_ellipse, sdf_half_plane, sdf_oriented_box, sdf_polygon, sdf_rounded_box, sdf_segment
from implicit_geometry_kernel.sdf.kernels.csg import blend, intersect, intersect_all, offset, shell, smooth_intersect, smooth_subtract, smooth_union, subtract, union, union_all
from implicit_geometry_kernel.sdf.kernels.transforms import SDFFrame, flapwise_aligned_web_angle, mirror_field_x, mirror_field_y, rotate_field, rotate_points, scale_field, translate_field
from implicit_geometry_kernel.sdf.kernels.grid import SDFGrid  # noqa: E402
from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF  # noqa: E402
from implicit_geometry_kernel.sdf.nodes import CallableField, Capsule, Circle, Box, Expr, HalfPlane, Intersect, Offset, OrientedBox, Polygon, RoundedBox, Rotate, Shell, Slab, Subtract, Union, compile_expr, eval_expr, eval_expr_batch, intersect_expr, subtract_expr, union_expr
from implicit_geometry_kernel.sdf.compile import compile_expr_to_node
from implicit_geometry_kernel.sdf.nodes.base import FeatureKind, FeatureSpec, FeatureTable, FieldSample, NO_FEATURE, SDFNode
from implicit_geometry_kernel.sdf.nodes.tagged import PriorityUnionNode, TagNode
from implicit_geometry_kernel.sdf.nodes.lofted import LoftedPolygonSDF
from implicit_geometry_kernel.sdf.nodes.span import ComplementNode, LawOffsetNode, MaxNode, SpanIntervalNode
from implicit_geometry_kernel.sdf.solve import GeometrySolver, SolveResult

__all__ = [
    "SDFCallable",
    # primitives
    "sdf_circle", "sdf_box", "sdf_half_plane", "sdf_segment", "sdf_capsule",
    "sdf_rounded_box", "sdf_ellipse", "sdf_oriented_box", "sdf_polygon",
    # csg
    "union", "intersect", "subtract", "offset", "shell", "blend",
    "smooth_union", "smooth_intersect", "smooth_subtract",
    "union_all", "intersect_all",
    # transforms
    "rotate_field", "translate_field", "scale_field",
    "mirror_field_x", "mirror_field_y",
    "rotate_points", "flapwise_aligned_web_angle", "SDFFrame",
    # grid
    "SDFGrid",
    # airfoil
    "AirfoilSDF",
    # ir
    "Expr", "Circle", "Box", "RoundedBox", "OrientedBox", "Capsule",
    "HalfPlane", "Polygon", "CallableField",
    "Union", "Intersect", "Subtract", "Offset", "Rotate", "Shell", "Slab",
    "eval_expr", "eval_expr_batch", "compile_expr", "compile_expr_to_node",
    "union_expr", "intersect_expr", "subtract_expr",
    # feature field
    "FieldSample", "FeatureTable", "FeatureSpec", "FeatureKind", "NO_FEATURE", "SDFNode",
    "TagNode", "PriorityUnionNode",
    "ComplementNode", "MaxNode", "LawOffsetNode", "SpanIntervalNode",
    "LoftedPolygonSDF", "GeometrySolver", "SolveResult",
]
