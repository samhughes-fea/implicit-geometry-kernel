"""Section geometry worker package public barrel.

G-row public products (see ``docs/reference/G_ROW_OWNERSHIP_MATRIX.md``):
``SectionGeometryRunner``, ``SectionGeometryResult``, ``MultiCellSection``,
``AirfoilSDF``, ``resolve_system_type`` / ``SYSTEM_TYPE_KEYS``, ``SectionLayout``.

Runner-adjacent helpers used by L4-EVAL stay in ``__all__``. Quadtree / hybrid
props backends are **G-PROPS debug** — preferred import is ``props_emit``;
they are intentionally **absent** from ``__all__`` (lazy shims remain for one
release). Plot / export helpers likewise demote to deep homes.
"""

from __future__ import annotations

from typing import Any

from .api import (
    SectionGeometryResult,
    SectionGeometryRunner,
    SectionPropertiesMeta,
    StationLayoutResolution,
    build_section_geometry,
    build_section_with_thicknesses,
    resolve_station_layout,
    resolve_system_type,
    section_geometry_prep_key,
)
from .constructive.labels import is_core_label
from .layout.section_layout import SectionLayout
from .layout.thickness_limits import (
    MIN_REALISTIC_SKIN_WALL_THICKNESS_M,
    MIN_REALISTIC_SPAR_WALL_THICKNESS_M,
    MIN_REALISTIC_WEB_WALL_THICKNESS_M,
    clamp_skin_thickness_m,
    clamp_spar_thickness_m,
    clamp_web_thickness_m,
)
from .layout.system_layout import SYSTEM_TYPE_KEYS, SystemLayoutSpec, layout_without_spar_caps
from .layout.spar_cap_runout import apply_spar_cap_runout, max_spar_from_cavity  # noqa: F401 — deep import
from .midlines.section_midline_host import SectionMidlineHost
from .sections import MultiCellSection
from .sdf.airfoil import AirfoilSDF

__version__ = "0.2.0"
__author__ = "Sam Hughes"

# G-row products + L4-EVAL runner-adjacent surface.
__all__ = [
    "SectionGeometryResult",
    "SectionGeometryRunner",
    "SectionPropertiesMeta",
    "StationLayoutResolution",
    "SectionLayout",
    "SectionMidlineHost",
    "build_section_geometry",
    "build_section_with_thicknesses",
    "section_geometry_prep_key",
    "AirfoilSDF",
    "MultiCellSection",
    "resolve_station_layout",
    "resolve_system_type",
    "SYSTEM_TYPE_KEYS",
    "SystemLayoutSpec",
    "layout_without_spar_caps",
    "is_core_label",
    "MIN_REALISTIC_SKIN_WALL_THICKNESS_M",
    "MIN_REALISTIC_SPAR_WALL_THICKNESS_M",
    "MIN_REALISTIC_WEB_WALL_THICKNESS_M",
    "clamp_skin_thickness_m",
    "clamp_spar_thickness_m",
    "clamp_web_thickness_m",
]

# One-release lazy shims (not in __all__). Prefer deep imports.
_LAZY: dict[str, tuple[str, str]] = {
    "apply_spar_cap_runout": (".layout.spar_cap_runout", "apply_spar_cap_runout"),
    "max_spar_from_cavity": (".layout.spar_cap_runout", "max_spar_from_cavity"),
    "build_section_view": (".api", "build_section_view"),
    "multicell_from_layout": (".sections", "multicell_from_layout"),
    "SDFCallable": (".layout.types", "SDFCallable"),
    "SDFGrid": (".sdf.grid", "SDFGrid"),
    "StructuralFamily": (".layout.structural", "StructuralFamily"),
    "FixedCapAnchor": (".layout.structural", "FixedCapAnchor"),
    "parse_structural_family": (".layout.structural", "parse_structural_family"),
    "parse_fixed_cap_anchor": (".layout.structural", "parse_fixed_cap_anchor"),
    "pitch_axis_x_from_le": (".layout.section_axes", "pitch_axis_x_from_le"),
    "max_thickness_chord_x": (".layout.section_axes", "max_thickness_chord_x"),
    "is_wall_solid_label": (".constructive.labels", "is_wall_solid_label"),
    "SectionPropertiesReport": (".props_emit.export", "SectionPropertiesReport"),
    "export_section_json": (".props_emit.export", "export_section_json"),
    "fiber_extents_from_airfoil": (".props_emit.props_handoff", "fiber_extents_from_airfoil"),
    "fiber_extents_from_outline": (".props_emit.props_handoff", "fiber_extents_from_outline"),
    "outer_polygon_from_airfoil": (".props_emit.props_handoff", "outer_polygon_from_airfoil"),
    "topology_for_section_properties": (".props_emit.props_handoff", "topology_for_section_properties"),
    "plot_section": (".visualisation", "plot_section"),
    "plot_sdf_field": (".visualisation", "plot_sdf_field"),
    # G-PROPS debug — prefer props_emit.*; kept lazy only for one release.
    "props_backend_from_env": (".props_emit.props_backend", "props_backend_from_env"),
    "quadtree_spec_for_airfoil": (".props_emit.props_backend", "quadtree_spec_for_airfoil"),
    "build_shared_section_quadtree": (".props_emit.area_mesh", "build_shared_section_quadtree"),
}

# Symbols banned from the public barrel forever (gate asserts absence from __all__).
_PROPS_DEBUG_BARREL_BAN = frozenset(
    {
        "props_backend_from_env",
        "quadtree_spec_for_airfoil",
        "build_shared_section_quadtree",
    }
)


def __getattr__(name: str) -> Any:
    if name in _LAZY:
        mod_name, attr = _LAZY[name]
        import importlib

        mod = importlib.import_module(mod_name, __name__)
        value = getattr(mod, attr)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__) | set(_LAZY))
