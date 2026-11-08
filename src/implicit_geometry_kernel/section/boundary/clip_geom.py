"""Clip polygons and capsules used by constructive emitters."""

from __future__ import annotations

import numpy as np
from implicit_geometry_kernel.section.boundary.polygons import (
    LineString,
    MultiPolygon,
    OPS,
    Polygon,
    airfoil_polygon,
    ensure_polygon,
    offset_polygon_inward,
)


def stadium_polygon(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    radius: float,
    *,
    resolution: int = 16,
) -> Polygon | MultiPolygon:
    """Capsule / stadium as a buffered line segment."""
    r = float(radius)
    if r <= 0:
        raise ValueError("stadium radius must be positive.")
    line = LineString([(float(x0), float(y0)), (float(x1), float(y1))])
    return ensure_polygon(
        OPS.buffer(line, r, quad_segs=resolution, cap_style=1, join_style=1)
    )


def half_plane_clip(geom, *, surface: str):
    """Clip to y>=0 (upper) or y<=0 (lower) in the chord frame."""
    minx, miny, maxx, maxy = geom.bounds
    pad = max(maxx - minx, maxy - miny, 1.0) + 1.0
    if surface == "upper":
        box = Polygon(
            [
                (minx - pad, 0.0),
                (maxx + pad, 0.0),
                (maxx + pad, maxy + pad),
                (minx - pad, maxy + pad),
            ]
        )
    elif surface == "lower":
        box = Polygon(
            [
                (minx - pad, miny - pad),
                (maxx + pad, miny - pad),
                (maxx + pad, 0.0),
                (minx - pad, 0.0),
            ]
        )
    else:
        raise ValueError(f"surface must be 'upper' or 'lower', got {surface!r}")
    return ensure_polygon(OPS.intersection(geom, box), allow_empty=True)


def chordwise_box(geom, x_start: float, x_end: float) -> Polygon | MultiPolygon:
    minx, miny, maxx, maxy = geom.bounds
    pad = max(maxy - miny, 1.0) + 1.0
    return Polygon(
        [
            (float(x_start), miny - pad),
            (float(x_end), miny - pad),
            (float(x_end), maxy + pad),
            (float(x_start), maxy + pad),
        ]
    )


def _padded_world_box(geom, pad_scale: float = 2.0) -> Polygon | MultiPolygon:
    minx, miny, maxx, maxy = geom.bounds
    span = max(maxx - minx, maxy - miny, 1.0)
    pad = pad_scale * span + 1.0
    return Polygon(
        [
            (minx - pad, miny - pad),
            (maxx + pad, miny - pad),
            (maxx + pad, maxy + pad),
            (minx - pad, maxy + pad),
        ]
    )


def flapwise_bay_clip(
    geom,
    *,
    nx: float,
    ny: float,
    clip_mode: str,
    c_lo: float | None,
    c_hi: float | None,
) -> Polygon | MultiPolygon:
    """Bay clip for flapwise webs: half_lt / half_gt / strip on n·p."""
    world = _padded_world_box(geom, pad_scale=4.0)
    minx, miny, maxx, maxy = world.bounds
    n = np.array([float(nx), float(ny)], dtype=float)
    n_norm = float(np.linalg.norm(n))
    if n_norm < 1e-15:
        raise ValueError("flapwise normal must be non-zero.")
    n = n / n_norm
    t = np.array([-n[1], n[0]], dtype=float)
    corners = np.array(
        [[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]], dtype=float
    )
    n_proj = corners @ n
    t_proj = corners @ t
    t0, t1 = float(np.min(t_proj)) - 1.0, float(np.max(t_proj)) + 1.0

    if clip_mode == "half_lt":
        if c_lo is None:
            raise ValueError("half_lt requires flapwise_c_lo")
        n0, n1 = float(np.min(n_proj)) - 1.0, float(c_lo)
    elif clip_mode == "half_gt":
        if c_lo is None:
            raise ValueError("half_gt requires flapwise_c_lo (bound)")
        n0, n1 = float(c_lo), float(np.max(n_proj)) + 1.0
    elif clip_mode == "strip":
        if c_lo is None or c_hi is None:
            raise ValueError("strip requires flapwise_c_lo and flapwise_c_hi")
        n0, n1 = float(min(c_lo, c_hi)), float(max(c_lo, c_hi))
    else:
        raise ValueError(f"Unknown flapwise_clip_mode {clip_mode!r}")

    nt_corners = [(n0, t0), (n1, t0), (n1, t1), (n0, t1)]
    xy = [tuple(n * nc + t * tc) for nc, tc in nt_corners]
    clip = Polygon(xy)
    return ensure_polygon(
        OPS.intersection(clip, OPS.envelope(minx, miny, maxx, maxy)), allow_empty=True
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
    "cavity_polygon",
    "chordwise_box",
    "flapwise_bay_clip",
    "half_plane_clip",
    "stadium_polygon",
]
