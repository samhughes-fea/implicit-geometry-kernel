"""Use-case-agnostic section geometry worker facade.

Mirrors the FE worker contract: no sizing/fatigue modes; callers choose
output gates (compute-cost knobs) and interpret results in orchestration.

Prep helpers live in ``prep``; assemble and emit live in ``build``.
"""

from implicit_geometry_kernel.section.spec.build import (
    _emit_section,
    build_section_geometry,
)
from implicit_geometry_kernel.section.spec.layout_resolution import (
    StationLayoutResolution,
    _as_section_layout,
    _layout_web_fracs,
    _reattach_catalogue_preset,
    _step_down_section_layout,
    resolve_station_layout,
    resolve_web_x_positions,
)
from implicit_geometry_kernel.section.spec.prep import (
    _compat_kwargs_from_zone_gauges,
    build_section_with_thicknesses,
    section_geometry_prep_key,
)
from implicit_geometry_kernel.section.spec.results import (
    ALL_OUTPUT_GATES,
    DroppedWeb,
    LayoutInput,
    OutputGate,
    SectionGeometryResult,
    SectionPropertiesMeta,
    normalize_output_gates,
)

__all__ = [
    "ALL_OUTPUT_GATES",
    "DroppedWeb",
    "LayoutInput",
    "OutputGate",
    "SectionGeometryResult",
    "SectionPropertiesMeta",
    "StationLayoutResolution",
    "_as_section_layout",
    "_compat_kwargs_from_zone_gauges",
    "_emit_section",
    "_layout_web_fracs",
    "_reattach_catalogue_preset",
    "_step_down_section_layout",
    "build_section_geometry",
    "build_section_with_thicknesses",
    "normalize_output_gates",
    "resolve_station_layout",
    "resolve_web_x_positions",
    "section_geometry_prep_key",
]
