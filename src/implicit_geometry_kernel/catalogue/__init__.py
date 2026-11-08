"""WX-Y-Z-NR catalogue: tokens, ``SystemLayoutSpec``, and the web-count cascade.

A token is a blade concept. This package holds the 2-D spec only (no stations).
Author a full blade through ``blade.concepts.BladeConcept``. Assemble one station
with ``SectionGeometry``.
"""

from implicit_geometry_kernel.catalogue.structural import (
    StructuralFamily,
    parse_structural_family,
)
from implicit_geometry_kernel.catalogue.system_layout import (
    SYSTEM_TYPE_KEYS,
    SYSTEM_TYPE_LAYOUT_KEYS,
    SystemLayoutSpec,
    WebAnchoring,
    WebOrientation,
    WebPositioning,
    effective_layout_for_web_count,
    family_min_webs,
    layout_without_spar_caps,
    resolve_system_type,
    step_down_cascade,
)

__all__ = [
    "SYSTEM_TYPE_KEYS",
    "SYSTEM_TYPE_LAYOUT_KEYS",
    "StructuralFamily",
    "SystemLayoutSpec",
    "WebAnchoring",
    "WebOrientation",
    "WebPositioning",
    "effective_layout_for_web_count",
    "family_min_webs",
    "layout_without_spar_caps",
    "parse_structural_family",
    "resolve_system_type",
    "step_down_cascade",
]
