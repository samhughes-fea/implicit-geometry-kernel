"""Midline-facing offset-ring helpers (re-export from package offset_math)."""

from implicit_geometry_kernel.section_geometry.offset_math import (
    OFFSET_JOIN_STYLE,
    OFFSET_MITRE_LIMIT,
    OFFSET_QUAD_SEGS,
    airfoil_offset_ring,
    locus_at_depth,
)

__all__ = [
    "OFFSET_JOIN_STYLE",
    "OFFSET_MITRE_LIMIT",
    "OFFSET_QUAD_SEGS",
    "airfoil_offset_ring",
    "locus_at_depth",
]
