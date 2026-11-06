"""Constructive wall-solid / core polygons and hybrid property helpers."""

from implicit_geometry_kernel.section.constructive.export import emit_constructive_polygons, emit_core_polygons, emit_wall_solid_polygons
from implicit_geometry_kernel.section.constructive.labels import INSERTS_CONSTRUCTIVE_BY_DEFAULT, classify_property_labels, is_constructive_label, is_core_label, is_insert_label, is_residual_label, is_skin_label, is_wall_solid_label
from implicit_geometry_kernel.section.boundary import OFFSET_JOIN_STYLE, OFFSET_MITRE_LIMIT, OFFSET_QUAD_SEGS, airfoil_offset_ring, airfoil_polygon, cavity_polygon, locus_at_depth, offset_polygon_inward
from implicit_geometry_kernel.section.constructive.properties import section_properties_from_polygons

__all__ = [
    "OFFSET_JOIN_STYLE",
    "OFFSET_MITRE_LIMIT",
    "OFFSET_QUAD_SEGS",
    "airfoil_offset_ring",
    "airfoil_polygon",
    "cavity_polygon",
    "classify_property_labels",
    "emit_constructive_polygons",
    "emit_core_polygons",
    "emit_wall_solid_polygons",
    "is_constructive_label",
    "INSERTS_CONSTRUCTIVE_BY_DEFAULT",
    "is_core_label",
    "is_insert_label",
    "is_residual_label",
    "is_skin_label",
    "is_wall_solid_label",
    "locus_at_depth",
    "offset_polygon_inward",
    "section_properties_from_polygons",
]
