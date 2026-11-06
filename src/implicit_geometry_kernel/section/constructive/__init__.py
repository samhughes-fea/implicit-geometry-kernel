"""Former home of the constructive polygons. Its modules moved to ``section.boundary``
(walls, cores, emitters, labels) and ``section.properties.polygons``; this package
keeps the old names as aliases until 0.6.0."""

from implicit_geometry_kernel.section.boundary.walls import emit_wall_solid_polygons
from implicit_geometry_kernel.section.boundary.cores import emit_constructive_polygons, emit_core_polygons
from implicit_geometry_kernel.section.boundary.labels import INSERTS_CONSTRUCTIVE_BY_DEFAULT, classify_property_labels, is_constructive_label, is_core_label, is_insert_label, is_residual_label, is_skin_label, is_wall_solid_label
from implicit_geometry_kernel.section.boundary import OFFSET_JOIN_STYLE, OFFSET_MITRE_LIMIT, OFFSET_QUAD_SEGS, airfoil_offset_ring, airfoil_polygon, cavity_polygon, locus_at_depth, offset_polygon_inward
from implicit_geometry_kernel.section.properties.polygons import section_properties_from_polygons

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
