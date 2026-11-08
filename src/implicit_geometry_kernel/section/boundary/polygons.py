"""Constructive polygon helpers; GEOS algorithms go through ``implicit_geometry_kernel.geos``.

Callers use ``PolygonOps`` and kernel geometry types (not Shapely).
Offsets use :func:`offset_polygon_inward` (mitre join; ``docs/section/bug_log.md``).
"""

from __future__ import annotations

import numpy as np

from implicit_geometry_kernel.contract.errors import InfeasibleSectionGeometry
from implicit_geometry_kernel.geos import lib as geos_lib
from implicit_geometry_kernel.geos.geom import LineString, MultiPolygon, Point, Polygon
from implicit_geometry_kernel.geos.wkb import dumps as _dumps
from implicit_geometry_kernel.geos.wkb import loads as _loads
from implicit_geometry_kernel.section.boundary.polygon_ops import (
    OFFSET_JOIN_STYLE,
    OFFSET_MITRE_LIMIT,
    OFFSET_QUAD_SEGS,
    PolygonOps,
    oriented_rings,
)


def _wkb(geom) -> bytes:
    wkb = getattr(geom, "wkb", None)
    if wkb is not None:
        return bytes(wkb)
    return _dumps(geom)


def _geom(wkb: bytes):
    return _loads(bytes(wkb))


def _map_xy(geom, fn):
    if geom is None or geom.is_empty:
        return geom
    if isinstance(geom, Polygon):
        holes = [fn(np.asarray(r.coords, float)) for r in geom.interiors]
        return Polygon(fn(np.asarray(geom.exterior.coords, float)), holes)
    if isinstance(geom, MultiPolygon):
        return MultiPolygon([_map_xy(g, fn) for g in geom.geoms])
    if isinstance(geom, LineString):
        return LineString(fn(np.asarray(geom.coords, float)))
    if isinstance(geom, Point):
        xy = fn(np.array([[float(geom.x), float(geom.y)]], float))[0]
        return Point(float(xy[0]), float(xy[1]))
    return geom


def _signed_area(xy: np.ndarray) -> float:
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


class GeosPolygonOps:
    """GEOS buffer and overlay via :mod:`implicit_geometry_kernel.geos.lib`."""

    def buffer(
        self,
        geom,
        distance: float,
        *,
        quad_segs: int | None = None,
        cap_style: int = 1,
        join_style: int | None = None,
        mitre_limit: float | None = None,
    ):
        return _geom(
            geos_lib.buffer(
                _wkb(geom),
                float(distance),
                quad_segs=OFFSET_QUAD_SEGS if quad_segs is None else int(quad_segs),
                cap_style=int(cap_style),
                join_style=OFFSET_JOIN_STYLE if join_style is None else int(join_style),
                mitre_limit=OFFSET_MITRE_LIMIT if mitre_limit is None else float(mitre_limit),
            )
        )

    def offset_inward(self, poly, depth: float, *, quad_segs: int | None = None):
        d = float(depth)
        if not np.isfinite(d) or d < 0.0:
            raise ValueError(f"offset_inward requires a finite depth >= 0; got {depth!r}.")
        return self.buffer(
            poly,
            -d,
            quad_segs=quad_segs,
            join_style=OFFSET_JOIN_STYLE,
            mitre_limit=OFFSET_MITRE_LIMIT,
        )

    def intersection(self, a, b):
        return _geom(geos_lib.intersection(_wkb(a), _wkb(b)))

    def difference(self, a, b):
        return _geom(geos_lib.difference(_wkb(a), _wkb(b)))

    def distance(self, a, b) -> float:
        return geos_lib.distance(_wkb(a), _wkb(b))

    def union(self, geoms):
        if geoms is None:
            return Polygon()
        if hasattr(geoms, "wkb") and not isinstance(geoms, (list, tuple)):
            return _geom(geos_lib.unary_union(_wkb(geoms)))
        seq = list(geoms)
        if not seq:
            return Polygon()
        return _geom(geos_lib.unary_union([_wkb(g) for g in seq]))

    def make_valid(self, geom):
        return _geom(geos_lib.make_valid(_wkb(geom)))

    def rotate(self, geom, angle_rad: float, *, cx: float = 0.0, cy: float = 0.0):
        if abs(float(angle_rad)) <= 1e-12:
            return geom
        c, s = np.cos(angle_rad), np.sin(angle_rad)
        ox, oy = float(cx), float(cy)

        def _rot(xy):
            x, y = xy[:, 0] - ox, xy[:, 1] - oy
            return np.column_stack([c * x - s * y + ox, s * x + c * y + oy])

        return _map_xy(geom, _rot)

    def envelope(self, minx, miny, maxx, maxy):
        return Polygon(
            [
                (float(minx), float(miny)),
                (float(maxx), float(miny)),
                (float(maxx), float(maxy)),
                (float(minx), float(maxy)),
            ]
        )

    def orient(self, poly, *, sign: float = 1.0):
        if not isinstance(poly, Polygon) or poly.is_empty:
            return poly
        want_ccw = float(sign) >= 0.0
        ext = np.asarray(poly.exterior.coords, float)
        if (_signed_area(ext) >= 0.0) != want_ccw:
            ext = ext[::-1]
        holes = []
        for ring in poly.interiors:
            xy = np.asarray(ring.coords, float)
            if (_signed_area(xy) >= 0.0) == want_ccw:
                xy = xy[::-1]
            holes.append(xy)
        return Polygon(ext, holes)

    def translate(self, geom, xoff: float, yoff: float):
        dx, dy = float(xoff), float(yoff)
        if dx == 0.0 and dy == 0.0:
            return geom
        return _map_xy(geom, lambda xy: xy + np.array([dx, dy], float))

    def explain_validity(self, geom) -> str:
        return geos_lib.is_valid_reason(_wkb(geom))


class KernelPolygonOps(GeosPolygonOps):
    """Numpy affine maps; offset and overlay stay GEOS."""


OPS = GeosPolygonOps()


def envelope(minx, miny, maxx, maxy):
    return OPS.envelope(minx, miny, maxx, maxy)


def offset_polygon_inward(poly, depth: float, *, quad_segs: int | None = None):
    """Inward normal offset of ``poly`` by ``depth``, the one offset operator."""
    d = float(depth)
    if not np.isfinite(d) or d < 0.0:
        raise ValueError(f"offset_polygon_inward requires a finite depth >= 0; got {depth!r}.")
    return OPS.offset_inward(poly, d, quad_segs=quad_segs)


def airfoil_polygon(airfoil_sdf) -> Polygon:
    verts = getattr(airfoil_sdf, "vertices", None)
    if verts is None:
        raise ValueError("AirfoilSDF-compatible object with vertices is required.")
    arr = np.asarray(verts, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2 or arr.shape[0] < 3:
        raise ValueError(f"Invalid airfoil vertices shape {arr.shape}.")
    poly = Polygon(arr)
    return validate_and_orient_polygon(poly, allow_empty=False)


def airfoil_offset_ring(
    airfoil_sdf, depth: float, *, quad_segs: int = OFFSET_QUAD_SEGS
) -> np.ndarray:
    """Vertices of the inward **normal**-offset locus ``phi_airfoil = -depth``."""
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


def locus_at_depth(
    airfoil_sdf, depth: float, *, quad_segs: int = OFFSET_QUAD_SEGS
) -> np.ndarray:
    return airfoil_offset_ring(airfoil_sdf, float(depth), quad_segs=quad_segs)


_MIN_PART_AREA = 1e-16


def validate_and_orient_polygon(poly, allow_empty: bool = False):
    """Validate and CCW-orient a polygon. Repair via buffer(0) if needed."""
    if allow_empty and (poly.is_empty or poly.area <= 0.0):
        return poly
    if not allow_empty and (poly.is_empty or poly.area <= 0.0):
        raise ValueError("Invalid polygon: empty or zero area")
    poly = OPS.orient(poly, sign=1.0)
    if not poly.is_valid:
        poly = OPS.buffer(poly, 0.0)
        if isinstance(poly, MultiPolygon):
            raise ValueError(
                "Polygon repair via buffer(0) produced MultiPolygon. "
                f"Original validity: {OPS.explain_validity(poly)}"
            )
        if allow_empty and poly.is_empty:
            return poly
    poly = OPS.orient(poly, sign=1.0)
    if not poly.is_valid or (not allow_empty and poly.is_empty):
        raise ValueError("Invalid polygon after orient")
    return poly


def unary_union(geoms):
    return OPS.union(geoms)


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
    try:
        if not geom.is_valid:
            geom = OPS.make_valid(geom)
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
        return ensure_polygon(OPS.union(parts), allow_empty=allow_empty)
    merged = OPS.union(geom)
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
    return OPS.rotate(geom, float(angle_rad), cx=float(cx), cy=float(cy))


from implicit_geometry_kernel.section.boundary.clip_geom import (  # noqa: E402
    cavity_polygon,
    chordwise_box,
    flapwise_bay_clip,
    half_plane_clip,
    stadium_polygon,
)

__all__ = [
    "LineString",
    "MultiPolygon",
    "OFFSET_JOIN_STYLE",
    "OFFSET_MITRE_LIMIT",
    "OFFSET_QUAD_SEGS",
    "PolygonOps",
    "Point",
    "Polygon",
    "GeosPolygonOps",
    "KernelPolygonOps",
    "OPS",
    "airfoil_offset_ring",
    "airfoil_polygon",
    "cavity_polygon",
    "chordwise_box",
    "ensure_polygon",
    "envelope",
    "flapwise_bay_clip",
    "half_plane_clip",
    "locus_at_depth",
    "offset_polygon_inward",
    "oriented_rings",
    "rotate_geom",
    "validate_and_orient_polygon",
    "stadium_polygon",
    "unary_union",
]
