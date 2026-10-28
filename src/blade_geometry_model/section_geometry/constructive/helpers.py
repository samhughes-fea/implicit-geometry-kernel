"""Shapely helpers for constructive wall-solid polygons.

Offset APIs (``offset_polygon_inward``, ``airfoil_offset_ring``,
``airfoil_polygon``) are **private implementation** of
:mod:`offset_locus`. Production code must import those names from
``offset_locus``, not from this module (enforced by
``tests/test_offset_import_boundary.py``).
"""

from __future__ import annotations

import numpy as np
from shapely import make_valid
from shapely.geometry import LineString, MultiPolygon, Polygon
from shapely.ops import unary_union

from blade_geometry_model.section_contract.errors import (
    InfeasibleSectionGeometry,
)
from blade_geometry_model.shoelace_integration import (
    validate_and_orient_polygon,
)

# Drop rings / parts below this absolute area (m^2). Thickened tip stations can
# produce Boolean slivers whose Shapely area is tiny but non-zero while shoelace
# reports A≈0 and previously crashed on Cx = Sx/(6A).
_MIN_PART_AREA = 1e-16

# The inward-offset join contract. Every offset locus in the package — skin
# mid-surface, inner cavity, cap slab, cap mid-surface — must come from these
# exact settings, because the package's junction exactness rests on loci at
# different depths being *parallel by construction*. A mitre join is the exact
# equidistant point for a locally-convex ring, so linear interpolation along an
# offset polyline stays on the locus; a different join style would not be
# equidistant, and two loci built with different settings would silently stop
# being parallel. Kept private: go through :func:`offset_polygon_inward`.
_OFFSET_QUAD_SEGS = 16
_OFFSET_JOIN_STYLE = 2  # mitre
_OFFSET_MITRE_LIMIT = 5.0


def offset_polygon_inward(poly, depth: float, *, quad_segs: int | None = None):
    """Inward normal offset of ``poly`` by ``depth``, the one offset operator.

    This is the single site that decides what "the offset locus at depth d"
    means. Prefer it over calling ``poly.buffer(-d, ...)``: the join settings
    are the package's parallelism contract (see ``docs/bug_log.md`` SG-001 and
    SG-003, where two loci that were not where each other expected produced
    10 mm displacements).

    Returns Shapely's result unchanged, **including an empty geometry** when
    the section is thinner than twice ``depth``. Callers decide whether empty
    is infeasible or simply means "no cavity here", and several rely on that.
    """
    d = float(depth)
    if not np.isfinite(d) or d < 0.0:
        raise ValueError(f"offset_polygon_inward requires a finite depth >= 0; got {depth!r}.")
    return poly.buffer(
        -d,
        quad_segs=_OFFSET_QUAD_SEGS if quad_segs is None else int(quad_segs),
        join_style=_OFFSET_JOIN_STYLE,
        mitre_limit=_OFFSET_MITRE_LIMIT,
    )


def airfoil_polygon(airfoil_sdf) -> Polygon:
    verts = getattr(airfoil_sdf, "vertices", None)
    if verts is None:
        raise ValueError("AirfoilSDF-compatible object with vertices is required.")
    arr = np.asarray(verts, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2 or arr.shape[0] < 3:
        raise ValueError(f"Invalid airfoil vertices shape {arr.shape}.")
    poly = Polygon(arr)
    return validate_and_orient_polygon(poly, allow_empty=False)


def airfoil_offset_ring(airfoil_sdf, depth: float, *, quad_segs: int = 16) -> np.ndarray:
    """Vertices of the inward **normal**-offset locus ``phi_airfoil = -depth``.

    Goes through :func:`offset_polygon_inward`, so this is the same offset
    family as :func:`emit_outer_skin_polygon` and :func:`emit_spar_cap_polygon`
    and a mid-surface built from this ring sits exactly inside the
    corresponding solid.

    Returns an ``(M, 2)`` ring **without** the repeated closing vertex.
    """
    d = float(depth)
    if not np.isfinite(d) or d <= 0.0:
        raise ValueError(f"airfoil_offset_ring requires depth > 0; got {depth!r}.")
    geom = offset_polygon_inward(airfoil_polygon(airfoil_sdf), d, quad_segs=quad_segs)
    if geom.is_empty:
        raise InfeasibleSectionGeometry(
            f"Inward offset at depth={d:g} is empty; the section is thinner than "
            "twice that depth."
        )
    if isinstance(geom, MultiPolygon):
        geom = max(geom.geoms, key=lambda g: g.area)
    if not isinstance(geom, Polygon) or geom.exterior is None:
        raise InfeasibleSectionGeometry(
            f"Inward offset at depth={d:g} did not yield a Polygon "
            f"(got {type(geom).__name__})."
        )
    ring = np.asarray(geom.exterior.coords, dtype=float)
    if ring.shape[0] >= 2 and float(np.hypot(*(ring[0] - ring[-1]))) <= 1e-15:
        ring = ring[:-1]
    if ring.shape[0] < 3:
        raise InfeasibleSectionGeometry(
            f"Inward offset at depth={d:g} produced too few vertices."
        )
    return ring


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


def chordwise_bay_clip(
    geom,
    *,
    x_start: float | None,
    x_end: float | None,
) -> Polygon:
    """Axis-aligned bay clip matching SDF HalfPlane chords (open ends allowed)."""
    box = _padded_world_box(geom)
    minx, miny, maxx, maxy = box.bounds
    x0 = float(x_start) if x_start is not None else minx
    x1 = float(x_end) if x_end is not None else maxx
    if x1 < x0:
        x0, x1 = x1, x0
    return Polygon(
        [
            (x0, miny),
            (x1, miny),
            (x1, maxy),
            (x0, maxy),
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
