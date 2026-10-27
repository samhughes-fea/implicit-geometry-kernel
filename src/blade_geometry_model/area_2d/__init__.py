"""2-D area meshing (quadtree default)."""

from .build import build_area_mesh_2d, build_area_meshes_2d
from .contour import ensure_closed, smooth_contour
from .quadtree import (
    QuadCell,
    QuadTreeMesh,
    build_root_quadtree,
    build_shared_quadtree,
    build_square_root_forest,
    contour_from_quadtree_marching_squares,
    refine_quadtree,
    refine_quadtree_multi,
    summarize_quadtree,
)
from .quadrature import (
    section_properties_all_labels,
    section_properties_from_quadtree,
)
from .spec import QuadTreeSpec, QuadtreeSpec
from .types import AreaMesh2D, FieldMesh2D, PhiSampler

__all__ = [
    "AreaMesh2D",
    "FieldMesh2D",
    "PhiSampler",
    "QuadCell",
    "QuadTreeMesh",
    "QuadTreeSpec",
    "QuadtreeSpec",
    "build_area_mesh_2d",
    "build_area_meshes_2d",
    "build_root_quadtree",
    "build_shared_quadtree",
    "build_square_root_forest",
    "contour_from_quadtree_marching_squares",
    "ensure_closed",
    "refine_quadtree",
    "refine_quadtree_multi",
    "section_properties_all_labels",
    "section_properties_from_quadtree",
    "smooth_contour",
    "summarize_quadtree",
]
