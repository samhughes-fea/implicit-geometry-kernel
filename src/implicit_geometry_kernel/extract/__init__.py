"""In-plane geometry extractors: polylines and feature curves, not meshes."""

from implicit_geometry_kernel.extract.curves import (
    classify_trailing_edge,
    feature_curve,
    locate_shell_vertex,
)
from implicit_geometry_kernel.extract.marching import extract_zero_loops

__all__ = [
    "classify_trailing_edge",
    "extract_zero_loops",
    "feature_curve",
    "locate_shell_vertex",
]
