"""Canonical offset-locus API for section_geometry.

Helpers in :mod:`helpers` own the Shapely implementation. This module is the
public facade: re-exports, documented join-contract constants, and thin
depth/cavity entry points used by clips, midlines, and emitters.
"""

from __future__ import annotations

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

from .helpers import (
    airfoil_offset_ring as _airfoil_offset_ring,
    airfoil_polygon as _airfoil_polygon,
    ensure_polygon,
    offset_polygon_inward as _offset_polygon_inward,
    rotate_geom,
)

# Public names for the inward-offset join contract (see helpers / bug_log SG-001).
OFFSET_QUAD_SEGS = 16
OFFSET_JOIN_STYLE = 2  # mitre
OFFSET_MITRE_LIMIT = 5.0


def offset_polygon_inward(poly, depth: float, *, quad_segs: int | None = None):
    """Inward normal offset; see :func:`helpers.offset_polygon_inward`."""
    return _offset_polygon_inward(poly, depth, quad_segs=quad_segs)


def airfoil_polygon(airfoil_sdf) -> Polygon:
    """Oriented airfoil contour polygon from SDF vertices."""
    return _airfoil_polygon(airfoil_sdf)


def airfoil_offset_ring(
    airfoil_sdf, depth: float, *, quad_segs: int = OFFSET_QUAD_SEGS
) -> np.ndarray:
    """Vertices of the inward normal-offset locus ``phi_airfoil = -depth``."""
    return _airfoil_offset_ring(airfoil_sdf, depth, quad_segs=quad_segs)


def locus_at_depth(
    airfoil_sdf, depth: float, *, quad_segs: int = OFFSET_QUAD_SEGS
) -> np.ndarray:
    """Ring vertices of the offset locus at ``depth`` (no repeated closing vertex)."""
    return airfoil_offset_ring(airfoil_sdf, float(depth), quad_segs=quad_segs)


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