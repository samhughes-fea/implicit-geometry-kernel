"""Area and moment integrals of 2-D regions.

Adaptive quadtree meshes over signed-distance fields with quadrature, and shoelace
(Green's theorem) formulas for polygons.
"""

from implicit_geometry_kernel.integration.build import build_area_mesh_2d, build_area_meshes_2d
from implicit_geometry_kernel.integration.contour import ensure_closed, smooth_contour
from implicit_geometry_kernel.integration.quadtree.tree import (
    QuadCell,
    QuadTreeMesh,
    build_root_quadtree,
    build_square_root_forest,
)
from implicit_geometry_kernel.integration.quadtree.soa import build_shared_quadtree
from implicit_geometry_kernel.integration.quadtree.marching import (
    contour_from_quadtree_marching_squares,
    summarize_quadtree,
)
from implicit_geometry_kernel.integration.quadtree.quadrature import (
    section_properties_all_labels,
    section_properties_from_quadtree,
)
from implicit_geometry_kernel.integration.shoelace import ring_sums_S, section_properties_S, section_properties_Sprime, validate_and_orient_polygon
from implicit_geometry_kernel.integration.spec import QuadTreeSpec, QuadtreeSpec
from implicit_geometry_kernel.integration.types import AreaMesh2D, FieldMesh2D, PhiSampler

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
    "ring_sums_S",
    "section_properties_S",
    "section_properties_Sprime",
    "validate_and_orient_polygon",
    "section_properties_all_labels",
    "section_properties_from_quadtree",
    "smooth_contour",
    "summarize_quadtree",
]
