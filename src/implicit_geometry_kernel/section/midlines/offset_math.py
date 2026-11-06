"""Moved to :mod:`implicit_geometry_kernel.section.boundary` (alias until 0.6.0)."""

from implicit_geometry_kernel.section.boundary import (  # noqa: F401
    OFFSET_JOIN_STYLE,
    OFFSET_MITRE_LIMIT,
    OFFSET_QUAD_SEGS,
    airfoil_offset_ring,
    airfoil_polygon,
    locus_at_depth,
    offset_polygon_inward,
)

__all__ = ["OFFSET_JOIN_STYLE", "OFFSET_MITRE_LIMIT", "OFFSET_QUAD_SEGS", "airfoil_offset_ring", "airfoil_polygon", "locus_at_depth", "offset_polygon_inward"]
