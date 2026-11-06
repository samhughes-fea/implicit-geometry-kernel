"""Former path of the quadtree quadrature. Moved to ``integration.quadtree.moments`` and
``integration.quadtree.quadrature``; this alias stays until 0.6.0."""

from implicit_geometry_kernel.integration.quadtree.moments import (
    _BOUNDARY_RECIPES,
    _accumulate_cell,
    _boundary_moments_batch,
    _boundary_recipe,
    _child_cells,
    _inside_polygon_from_cell,
    _leaf_corner_coords,
    _leaf_edge_crossings,
    _polygon_moments,
    _rectangle_moments,
    _rectangle_moments_batch,
    _shoelace_moments_batch,
)
from implicit_geometry_kernel.integration.quadtree.quadrature import (
    _accumulate_moments_vectorized,
    _active_leaf_indices,
    _attach_center_phi,
    _label_corner_phi,
    _raw_to_section_properties,
    label_boundary_segments,
    label_fill_polygons,
    section_properties_all_labels,
    section_properties_from_quadtree,
)

__all__ = [
    "_BOUNDARY_RECIPES",
    "_accumulate_cell",
    "_boundary_moments_batch",
    "_boundary_recipe",
    "_child_cells",
    "_inside_polygon_from_cell",
    "_leaf_corner_coords",
    "_leaf_edge_crossings",
    "_polygon_moments",
    "_rectangle_moments",
    "_rectangle_moments_batch",
    "_shoelace_moments_batch",
    "_accumulate_moments_vectorized",
    "_active_leaf_indices",
    "_attach_center_phi",
    "_label_corner_phi",
    "_raw_to_section_properties",
    "label_boundary_segments",
    "label_fill_polygons",
    "section_properties_all_labels",
    "section_properties_from_quadtree",
]
