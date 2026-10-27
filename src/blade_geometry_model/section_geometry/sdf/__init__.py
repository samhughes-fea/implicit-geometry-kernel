"""
blade_geometry_model.section_geometry.sdf
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
    AirfoilSDF

NACA generators:
    naca_four_digit_vertices, naca_five_digit_vertices, naca_six_series_vertices,
    airfoil_vertices_from_spanwise, spanwise_airfoil_label

CSG IR:
    Expr, Circle, Box, RoundedBox, OrientedBox, Capsule, HalfPlane, Polygon,
    CallableField, Union, Intersect, Subtract, Offset, Rotate,
    eval_expr, compile_expr

Cache:
    SectionEvalCache

Protocol:
    SDFCallable
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable

@runtime_checkable
class SDFCallable(Protocol):
    """Callable signed-distance field interface: (x, y) -> ndarray."""

    def __call__(self, x: Any, y: Any) -> Any: ...

from .primitives import (  # noqa: E402
    sdf_box,
    sdf_capsule,
    sdf_circle,
    sdf_ellipse,
    sdf_half_plane,
    sdf_oriented_box,
    sdf_polygon,
    sdf_rounded_box,
    sdf_segment,
)
from .csg import (  # noqa: E402
    blend,
    intersect,
    intersect_all,
    offset,
    shell,
    smooth_intersect,
    smooth_subtract,
    smooth_union,
    subtract,
    union,
    union_all,
)
from .transforms import (  # noqa: E402
    SDFFrame,
    flapwise_aligned_web_angle,
    mirror_field_x,
    mirror_field_y,
    rotate_field,
    rotate_points,
    scale_field,
    translate_field,
)
from .grid import SDFGrid  # noqa: E402
from .airfoil import AirfoilSDF  # noqa: E402
from .naca import (  # noqa: E402
    airfoil_vertices_from_spanwise,
    naca_five_digit_vertices,
    naca_four_digit_vertices,
    naca_six_series_vertices,
    spanwise_airfoil_label,
)
from .ir import (  # noqa: E402
    CallableField,
    Capsule,
    Circle,
    Box,
    Expr,
    HalfPlane,
    Intersect,
    Offset,
    OrientedBox,
    Polygon,
    RoundedBox,
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
from .cache import SectionEvalCache  # noqa: E402

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
    # naca
    "naca_four_digit_vertices", "naca_five_digit_vertices",
    "naca_six_series_vertices", "airfoil_vertices_from_spanwise",
    "spanwise_airfoil_label",
    # ir
    "Expr", "Circle", "Box", "RoundedBox", "OrientedBox", "Capsule",
    "HalfPlane", "Polygon", "CallableField",
    "Union", "Intersect", "Subtract", "Offset", "Rotate", "Shell", "Slab",
    "eval_expr", "eval_expr_batch", "compile_expr",
    "union_expr", "intersect_expr", "subtract_expr",
    # cache
    "SectionEvalCache",
]
