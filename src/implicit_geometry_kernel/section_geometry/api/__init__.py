"""Public API for section geometry."""

from .factory import (
    SYSTEM_TYPE_KEYS,
    build_section_view,
    resolve_system_type,
)
from .runner import SectionGeometryRunner
from .worker import (
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
from ..props_emit.props_backend import (
    props_backend_from_env,
    quadtree_spec_for_airfoil,
    quadtree_spec_to_dict,
)

__all__ = [
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
