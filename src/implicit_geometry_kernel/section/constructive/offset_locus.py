"""G-POLY offset-locus facade: re-exports G-MID offset math plus cavity helpers."""

from __future__ import annotations

from shapely.geometry import MultiPolygon, Polygon

from implicit_geometry_kernel.section.constructive.helpers import ensure_polygon, rotate_geom
from implicit_geometry_kernel.section.midlines.offset_math import (
    OFFSET_JOIN_STYLE,
    OFFSET_MITRE_LIMIT,
    OFFSET_QUAD_SEGS,
    airfoil_offset_ring,
    airfoil_polygon,
    locus_at_depth,
    offset_polygon_inward,
)


def cavity_polygon(
    airfoil_sdf,
    skin_thickness: float,
    *,
    airfoil_poly: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon:
    """Inner-skin cavity: Offset(airfoil, -skin_thickness)."""
    poly = airfoil_poly if airfoil_poly is not None else airfoil_polygon(airfoil_sdf)
    cavity = offset_polygon_inward(poly, float(skin_thickness))
    return ensure_polygon(cavity, allow_empty=True)


__all__ = [
    "OFFSET_JOIN_STYLE",
    "OFFSET_MITRE_LIMIT",
    "OFFSET_QUAD_SEGS",
    "airfoil_offset_ring",
    "airfoil_polygon",
    "cavity_polygon",
    "ensure_polygon",
    "locus_at_depth",
    "offset_polygon_inward",
    "rotate_geom",
]
