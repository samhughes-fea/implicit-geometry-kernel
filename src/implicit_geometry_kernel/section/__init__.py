"""Section geometry worker package public barrel.

Public products (see ``docs/section/reference/G_ROW_OWNERSHIP_MATRIX.md``):
``SectionGeometry``, ``SectionGeometryResult``, ``MultiCellSection``,
``resolve_system_type`` / ``SYSTEM_TYPE_KEYS``, ``SectionLayout``.

Prefer ``Airfoil.to_sdf()`` (or ``sdf.shapes.airfoil.AirfoilSDF``) for implicit
airfoils. Quadtree / hybrid props backends are debug — import from
``section.properties``, not this barrel. Thickness clamps live in
``section.layout.thickness_limits``. ``build_section_geometry`` is
``section.spec``.
"""

from __future__ import annotations

from implicit_geometry_kernel.section.spec import (
    Gauges,
    SectionGeometry,
    SectionGeometryResult,
    SectionOptions,
    SectionPropertiesMeta,
    SectionSpec,
    StationLayoutResolution,
    resolve_station_layout,
    resolve_system_type,
)
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
from implicit_geometry_kernel.catalogue import (
    SYSTEM_TYPE_KEYS,
    SYSTEM_TYPE_LAYOUT_KEYS,
    SystemLayoutSpec,
)
from implicit_geometry_kernel.section.assembly import MultiCellSection

__version__ = "0.8.2"
__author__ = "Sam Hughes"

__all__ = [
    "Gauges",
    "SectionOptions",
    "SectionSpec",
    "SectionGeometry",
    "SectionGeometryResult",
    "SectionPropertiesMeta",
    "StationLayoutResolution",
    "SectionLayout",
    "MultiCellSection",
    "resolve_station_layout",
    "resolve_system_type",
    "SYSTEM_TYPE_KEYS",
    "SYSTEM_TYPE_LAYOUT_KEYS",
    "SystemLayoutSpec",
]
