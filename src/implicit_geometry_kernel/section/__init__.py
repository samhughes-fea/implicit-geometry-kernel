"""Section geometry worker package public barrel.

G-row public products (see ``docs/section/reference/G_ROW_OWNERSHIP_MATRIX.md``):
``SectionGeometry``, ``SectionGeometryResult``, ``MultiCellSection``,
``resolve_system_type`` / ``SYSTEM_TYPE_KEYS``, ``SectionLayout``.

Prefer ``Airfoil.to_sdf()`` (or ``sdf.airfoil.AirfoilSDF``) for implicit
airfoils. Quadtree / hybrid props backends are G-PROPS debug — import from
``section.properties``, not this barrel. Thickness clamps live in
``section.layout.thickness_limits``. ``build_section_geometry`` is
``section.api``.
"""

from __future__ import annotations

from implicit_geometry_kernel.section.api import (
    Gauges,
    SectionGeometry,
    SectionGeometryResult,
    SectionGeometryRunner,
    SectionOptions,
    SectionPropertiesMeta,
    SectionSpec,
    StationLayoutResolution,
    resolve_station_layout,
    resolve_system_type,
)
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
from implicit_geometry_kernel.section.layout.system_layout import SYSTEM_TYPE_KEYS, SystemLayoutSpec
from implicit_geometry_kernel.section.assembly import MultiCellSection

__version__ = "0.4.0"
__author__ = "Sam Hughes"

__all__ = [
    "Gauges",
    "SectionOptions",
    "SectionSpec",
    "SectionGeometry",
    "SectionGeometryResult",
    "SectionGeometryRunner",
    "SectionPropertiesMeta",
    "StationLayoutResolution",
    "SectionLayout",
    "MultiCellSection",
    "resolve_station_layout",
    "resolve_system_type",
    "SYSTEM_TYPE_KEYS",
    "SystemLayoutSpec",
]
