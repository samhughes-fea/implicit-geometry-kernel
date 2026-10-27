"""Shared offset-locus math for constructive and midlines.

Neutral package home so G-MID does not import constructive.
"""

from __future__ import annotations

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

from blade_geometry_model.section_geometry.constructive.helpers import (
    airfoil_offset_ring as _airfoil_offset_ring,
    airfoil_polygon as _airfoil_polygon,
    ensure_polygon,
    offset_polygon_inward as _offset_polygon_inward,
)

OFFSET_QUAD_SEGS = 16
OFFSET_JOIN_STYLE = 2  # mitre
OFFSET_MITRE_LIMIT = 5.0


def offset_polygon_inward(poly, depth: float, *, quad_segs: int | None = None):
    return _offset_polygon_inward(poly, depth, quad_segs=quad_segs)


def airfoil_polygon(airfoil_sdf) -> Polygon:
    return _airfoil_polygon(airfoil_sdf)


def airfoil_offset_ring(
    airfoil_sdf, depth: float, *, quad_segs: int = OFFSET_QUAD_SEGS
) -> np.ndarray:
    return _airfoil_offset_ring(airfoil_sdf, depth, quad_segs=quad_segs)


def locus_at_depth(
    airfoil_sdf, depth: float, *, quad_segs: int = OFFSET_QUAD_SEGS
) -> np.ndarray:
    return airfoil_offset_ring(airfoil_sdf, float(depth), quad_segs=quad_segs)


def cavity_polygon(
    airfoil_sdf,
    skin_thickness: float,
    *,
    airfoil_poly: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon:
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
    "locus_at_depth",
    "offset_polygon_inward",
]
