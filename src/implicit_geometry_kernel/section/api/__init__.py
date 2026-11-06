"""Public API for section geometry."""

from implicit_geometry_kernel.section.api.factory import SYSTEM_TYPE_KEYS, build_section_view, resolve_system_type
from implicit_geometry_kernel.section.api.runner import SectionGeometry, SectionGeometryRunner
from implicit_geometry_kernel.section.api.worker import ALL_OUTPUT_GATES, DroppedWeb, SectionGeometryResult, SectionPropertiesMeta, StationLayoutResolution, build_section_geometry, build_section_with_thicknesses, normalize_output_gates, resolve_station_layout, section_geometry_prep_key
from implicit_geometry_kernel.section.properties.backend import props_backend_from_env, quadtree_spec_for_airfoil, quadtree_spec_to_dict

__all__ = [
    "SectionGeometry",
    "SectionGeometryRunner",
    "SYSTEM_TYPE_KEYS",
    "build_section_view",
    "resolve_system_type",
    "ALL_OUTPUT_GATES",
    "DroppedWeb",
    "SectionGeometryResult",
    "SectionPropertiesMeta",
    "StationLayoutResolution",
    "build_section_geometry",
    "build_section_with_thicknesses",
    "normalize_output_gates",
    "resolve_station_layout",
    "section_geometry_prep_key",
    "props_backend_from_env",
    "quadtree_spec_for_airfoil",
    "quadtree_spec_to_dict",
]
