"""Section input types and the geometry worker (validation, resolve, build)."""

from implicit_geometry_kernel.section.properties.backend import (
    props_backend_from_env,
    quadtree_spec_for_airfoil,
    quadtree_spec_to_dict,
)
from implicit_geometry_kernel.section.spec.factory import (
    SYSTEM_TYPE_KEYS,
    build_section_view,
    resolve_system_type,
)
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


def __getattr__(name: str):
    if name == "SectionGeometryRunner":
        import warnings

        warnings.warn(
            "SectionGeometryRunner is deprecated; use SectionGeometry. Removed in 0.6.0.",
            DeprecationWarning,
            stacklevel=2,
        )
        return SectionGeometry
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
