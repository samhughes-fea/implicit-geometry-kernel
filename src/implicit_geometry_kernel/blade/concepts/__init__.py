"""Blade concepts: catalogue token plus spanwise scheme, and web anchoring."""

from implicit_geometry_kernel.blade.concepts.concept import (
    BladeConcept,
    SpanwiseScheme,
    catalogue_ribs,
)
from implicit_geometry_kernel.blade.concepts.web_anchoring import (
    PaRegions,
    blade_x_le_te,
    lc_web_fractions,
    le_web_offsets_m,
    pa_regions,
    pa_web_x_from_pitch,
)
from implicit_geometry_kernel.catalogue import (
    SYSTEM_TYPE_KEYS,
    SYSTEM_TYPE_LAYOUT_KEYS,
    SystemLayoutSpec,
    resolve_system_type,
)

__all__ = [
    "BladeConcept",
    "PaRegions",
    "SYSTEM_TYPE_KEYS",
    "SYSTEM_TYPE_LAYOUT_KEYS",
    "SpanwiseScheme",
    "SystemLayoutSpec",
    "blade_x_le_te",
    "catalogue_ribs",
    "lc_web_fractions",
    "le_web_offsets_m",
    "pa_regions",
    "pa_web_x_from_pitch",
    "resolve_system_type",
]
