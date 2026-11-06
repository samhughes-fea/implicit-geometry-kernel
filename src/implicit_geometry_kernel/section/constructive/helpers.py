"""Shapely helpers for constructive wall-solid polygons.

Offset APIs live in ``section.midlines.offset_math``. This module re-exports
them for local constructive helpers. Production code must import those names
from ``offset_locus``, not from this module (enforced by
``tests/test_offset_import_boundary.py``).
"""

from __future__ import annotations

import numpy as np
from shapely import make_valid
from shapely.geometry import LineString, MultiPolygon, Polygon
from shapely.ops import unary_union

from implicit_geometry_kernel.integration.shoelace import (
    validate_and_orient_polygon,
)
from implicit_geometry_kernel.section.midlines.offset_math import (
    airfoil_offset_ring,
    airfoil_polygon,
    offset_polygon_inward,
)

# Drop rings / parts below this absolute area (m^2). Thickened tip stations can
# produce Boolean slivers whose Shapely area is tiny but non-zero while shoelace
# reports A≈0 and previously crashed on Cx = Sx/(6A).
_MIN_PART_AREA = 1e-16


def _positive_area(geom) -> bool:
    try:
        return float(getattr(geom, "area", 0.0) or 0.0) > _MIN_PART_AREA
    except Exception:
        return False


def ensure_polygon(geom, *, allow_empty: bool = False) -> Polygon | MultiPolygon:
    if geom is None or geom.is_empty:
        if allow_empty:
            return Polygon()
        raise ValueError("Empty geometry after constructive Boolean.")
    # Repair self-touches / bow-ties from offset+Boolean before orientation.
    try:
        if not geom.is_valid:
            geom = make_valid(geom)
    except Exception:
        pass
    if geom is None or geom.is_empty or not _positive_area(geom):
        if allow_empty:
            return Polygon()
        raise ValueError("Empty geometry after constructive Boolean.")
    if isinstance(geom, Polygon):
        return validate_and_orient_polygon(geom, allow_empty=allow_empty)
    if isinstance(geom, MultiPolygon):
        parts = []
        for g in geom.geoms:
            if g.is_empty or not _positive_area(g):
                continue
            try:
                parts.append(validate_and_orient_polygon(g, allow_empty=False))
            except ValueError:
                continue
        if not parts:
            if allow_empty:
                return Polygon()
            raise ValueError("MultiPolygon had no positive-area parts.")
        if len(parts) == 1:
            return parts[0]
        return MultiPolygon(parts)
    # GeometryCollection / LineString / etc. — extract polygonal parts once.
    geoms = getattr(geom, "geoms", None)
    if geoms is not None:
        parts = []
        for g in geoms:
            if isinstance(g, (Polygon, MultiPolygon)) and not g.is_empty and _positive_area(g):
                parts.append(g)
        if not parts:
            if allow_empty:
                return Polygon()
            raise ValueError("No polygonal parts in geometry collection.")
        return ensure_polygon(unary_union(parts), allow_empty=allow_empty)
    merged = unary_union(geom)
    if merged is geom or (
        hasattr(merged, "geom_type")
        and merged.geom_type not in {"Polygon", "MultiPolygon"}
    ):
        if allow_empty:
            return Polygon()
        raise ValueError(
            f"Cannot coerce geometry type {getattr(geom, 'geom_type', type(geom))} to polygon."
        )
    return ensure_polygon(merged, allow_empty=allow_empty)


def rotate_geom(geom, angle_rad: float, *, cx: float = 0.0, cy: float = 0.0):
    """Rotate geometry about (cx, cy) by ``angle_rad`` (radians, CCW)."""
    if abs(float(angle_rad)) <= 1e-12:
        return geom
    from shapely.affinity import rotate

    return rotate(geom, float(np.degrees(angle_rad)), origin=(float(cx), float(cy)))


def stadium_polygon(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    radius: float,
    *,
    resolution: int = 16,
) -> Polygon:
    """Capsule / stadium as a buffered line segment."""
    r = float(radius)
    if r <= 0:
        raise ValueError("stadium radius must be positive.")
    line = LineString([(float(x0), float(y0)), (float(x1), float(y1))])
    return ensure_polygon(line.buffer(r, quad_segs=resolution, cap_style=1, join_style=1))


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
    return ensure_polygon(geom.intersection(box), allow_empty=True)


def chordwise_box(geom, x_start: float, x_end: float) -> Polygon:
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


def _padded_world_box(geom, pad_scale: float = 2.0) -> Polygon:
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
) -> Polygon:
    """Bay clip for flapwise webs: half_lt / half_gt / strip on n·p.

    Matches ``_parallel_web_half_lt`` (n·p < c), ``_parallel_web_half_gt``
    (n·p > c), and ``_parallel_web_strip_clip`` (lo ≤ n·p ≤ hi).
    """
    from shapely.geometry import box as shapely_box

    world = _padded_world_box(geom, pad_scale=4.0)
    minx, miny, maxx, maxy = world.bounds
    # Dense sample of candidate clip polygon via rotated rectangle approximation:
    # build a large axis-aligned box in (n, t) coordinates then rotate.
    n = np.array([float(nx), float(ny)], dtype=float)
    n_norm = float(np.linalg.norm(n))
    if n_norm < 1e-15:
        raise ValueError("flapwise normal must be non-zero.")
    n = n / n_norm
    t = np.array([-n[1], n[0]], dtype=float)
    # World extent along n and t
    corners = np.array(
        [[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]], dtype=float
    )
    n_proj = corners @ n
    t_proj = corners @ t
    t0, t1 = float(t_proj.min()) - 1.0, float(t_proj.max()) + 1.0

    if clip_mode == "half_lt":
        if c_lo is None:
            raise ValueError("half_lt requires flapwise_c_lo")
        n0, n1 = float(n_proj.min()) - 1.0, float(c_lo)
    elif clip_mode == "half_gt":
        if c_lo is None:
            raise ValueError("half_gt requires flapwise_c_lo (bound)")
        n0, n1 = float(c_lo), float(n_proj.max()) + 1.0
    elif clip_mode == "strip":
        if c_lo is None or c_hi is None:
            raise ValueError("strip requires flapwise_c_lo and flapwise_c_hi")
        n0, n1 = float(min(c_lo, c_hi)), float(max(c_lo, c_hi))
    else:
        raise ValueError(f"Unknown flapwise_clip_mode {clip_mode!r}")

    # Rectangle corners in (n, t) then map back to xy: p = n*n_coord + t*t_coord
    nt_corners = [
        (n0, t0),
        (n1, t0),
        (n1, t1),
        (n0, t1),
    ]
    xy = [tuple(n * nc + t * tc) for nc, tc in nt_corners]
    clip = Polygon(xy)
    # Keep finite intersection with world (numerical safety)
    return ensure_polygon(clip.intersection(shapely_box(minx, miny, maxx, maxy)), allow_empty=True)


__all__ = [
    "airfoil_offset_ring",
    "airfoil_polygon",
    "chordwise_box",
    "ensure_polygon",
    "flapwise_bay_clip",
    "half_plane_clip",
    "offset_polygon_inward",
    "rotate_geom",
    "stadium_polygon",
]
