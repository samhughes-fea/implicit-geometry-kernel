"""Section input types and the geometry worker (validation, resolve, build)."""

from implicit_geometry_kernel.section.properties.backend import (
    props_backend_from_env,
    quadtree_spec_for_airfoil,
    quadtree_spec_to_dict,
)
from implicit_geometry_kernel.catalogue import SYSTEM_TYPE_KEYS, resolve_system_type
from implicit_geometry_kernel.section.assembly import build_section_view
from implicit_geometry_kernel.section.spec.options import SectionOptions
from implicit_geometry_kernel.section.spec.runner import (
    SectionGeometry,
)
from implicit_geometry_kernel.section.spec.section_spec import Gauges, SectionSpec
from implicit_geometry_kernel.section.spec.worker import (
    ALL_OUTPUT_GATES,
    DroppedWeb,
    SectionGeometryResult,
    SectionPropertiesMeta,
    StationLayoutResolution,
    build_section_geometry,
    build_section_with_thicknesses,
    normalize_output_gates,
    resolve_station_layout,
    section_geometry_prep_key,
)

__all__ = [
    "Gauges",
    "SectionOptions",
    "SectionSpec",
    "SectionGeometry",
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
