"""Moved to :mod:`implicit_geometry_kernel.section.boundary` (alias until 0.6.0)."""

from implicit_geometry_kernel.section.boundary import (  # noqa: F401
    airfoil_offset_ring,
    airfoil_polygon,
    chordwise_box,
    ensure_polygon,
    flapwise_bay_clip,
    half_plane_clip,
    offset_polygon_inward,
    rotate_geom,
    stadium_polygon,
)

__all__ = ["airfoil_offset_ring", "airfoil_polygon", "chordwise_box", "ensure_polygon", "flapwise_bay_clip", "half_plane_clip", "offset_polygon_inward", "rotate_geom", "stadium_polygon"]
