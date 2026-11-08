"""Numpy-backed 2-D geometries. Overlay lives on PolygonOps, not on these types."""

from __future__ import annotations

import numpy as np


def _as_xy(coords) -> np.ndarray:
    if coords is None:
        return np.zeros((0, 2), float)
    if isinstance(coords, CoordinateSequence):
        return coords._xy.copy()
    arr = np.asarray(coords, float)
    if arr.size == 0:
        return np.zeros((0, 2), float)
    if arr.ndim == 1 and arr.shape[0] == 2:
        return arr.reshape(1, 2)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(f"Expected (n, 2) coordinates, got {arr.shape}")
    return np.ascontiguousarray(arr[:, :2], dtype=float)


def _signed_area(xy: np.ndarray) -> float:
    if xy.shape[0] < 3:
        return 0.0
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _hausdorff(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0 or b.size == 0:
        return 0.0

    def _h(p, q):
        d = np.sqrt(((p[:, None, :] - q[None, :, :]) ** 2).sum(axis=2))
        return float(np.max(np.min(d, axis=1)))

    return max(_h(a, b), _h(b, a))


class CoordinateSequence:
    def __init__(self, xy):
        self._xy = _as_xy(xy)

    def __array__(self, dtype=None):
        return np.asarray(self._xy, dtype=dtype)

    def __len__(self):
        return int(self._xy.shape[0])

    def __iter__(self):
        for row in self._xy:
            yield (float(row[0]), float(row[1]))

    def __getitem__(self, i):
        row = self._xy[i]
        if isinstance(i, slice):
            return CoordinateSequence(row)
        return (float(row[0]), float(row[1]))

    @property
    def xy(self):
        return self._xy[:, 0].copy(), self._xy[:, 1].copy()

    def hausdorff_distance(self, other) -> float:
        b = other._xy if isinstance(other, CoordinateSequence) else _as_xy(getattr(other, "coords", other))
        return _hausdorff(self._xy, b)


class LinearRing:
    def __init__(self, coords):
        self.coords = CoordinateSequence(coords)

    @property
    def xy(self):
        return self.coords.xy

    def hausdorff_distance(self, other) -> float:
        o = other.coords if isinstance(other, LinearRing) else other
        return self.coords.hausdorff_distance(o)


class Point:
    geom_type = "Point"

    def __init__(self, x, y=None):
        if y is None:
            if isinstance(x, Point):
                self.x, self.y = float(x.x), float(x.y)
            else:
                xy = np.asarray(x, float).reshape(-1)
                self.x, self.y = float(xy[0]), float(xy[1])
        else:
            self.x, self.y = float(x), float(y)

    @property
    def xy(self):
        return (np.array([self.x]), np.array([self.y]))

    @property
    def is_empty(self) -> bool:
        return not (np.isfinite(self.x) and np.isfinite(self.y))

    @property
    def coords(self):
        return CoordinateSequence([[self.x, self.y]])

    @property
    def bounds(self):
        return (self.x, self.y, self.x, self.y)

    @property
    def wkb(self) -> bytes:
        from implicit_geometry_kernel.geos.wkb import dumps

        return dumps(self)


class LineString:
    geom_type = "LineString"

    def __init__(self, coords=None):
        self.coords = CoordinateSequence(coords)

    @property
    def is_empty(self) -> bool:
        return len(self.coords) == 0

    @property
    def bounds(self):
        xy = self.coords._xy
        if xy.size == 0:
            return (0.0, 0.0, 0.0, 0.0)
        return (float(xy[:, 0].min()), float(xy[:, 1].min()), float(xy[:, 0].max()), float(xy[:, 1].max()))

    @property
    def wkb(self) -> bytes:
        from implicit_geometry_kernel.geos.wkb import dumps

        return dumps(self)

    def distance(self, other) -> float:
        a = self.coords._xy
        if isinstance(other, Point):
            b = np.array([[other.x, other.y]], float)
        else:
            b = _as_xy(getattr(other, "coords", other))
        if a.size == 0 or b.size == 0:
            return 0.0
        d = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2))
        return float(np.min(d))


class Polygon:
    geom_type = "Polygon"

    def __init__(self, shell=None, holes=None):
        if isinstance(shell, Polygon):
            self._ext = LinearRing(shell.exterior.coords)
            self._holes = tuple(LinearRing(r.coords) for r in shell.interiors)
            return
        if shell is None:
            self._ext = LinearRing(np.zeros((0, 2), float))
            self._holes = ()
            return
        self._ext = LinearRing(shell)
        self._holes = tuple(LinearRing(h) for h in (holes or ()))

    @property
    def exterior(self) -> LinearRing:
        return self._ext

    @property
    def interiors(self):
        return self._holes

    @property
    def is_empty(self) -> bool:
        return len(self._ext.coords) < 3

    @property
    def area(self) -> float:
        a = abs(_signed_area(self._ext.coords._xy))
        for h in self._holes:
            a -= abs(_signed_area(h.coords._xy))
        return float(max(a, 0.0))

    @property
    def bounds(self):
        xy = self._ext.coords._xy
        if xy.size == 0:
            return (0.0, 0.0, 0.0, 0.0)
        return (float(xy[:, 0].min()), float(xy[:, 1].min()), float(xy[:, 0].max()), float(xy[:, 1].max()))

    @property
    def is_valid(self) -> bool:
        from implicit_geometry_kernel.geos import lib
        from implicit_geometry_kernel.geos.wkb import dumps

        if self.is_empty:
            return True
        return lib.is_valid(dumps(self))

    @property
    def centroid(self) -> Point:
        xy = self._ext.coords._xy
        if xy.shape[0] < 3:
            return Point(0.0, 0.0)
        a2 = _signed_area(xy) * 2.0
        if abs(a2) < 1e-18:
            return Point(float(xy[:, 0].mean()), float(xy[:, 1].mean()))
        x, y = xy[:, 0], xy[:, 1]
        cx = np.sum((x + np.roll(x, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)) / (3.0 * a2)
        cy = np.sum((y + np.roll(y, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)) / (3.0 * a2)
        return Point(float(cx), float(cy))

    def representative_point(self) -> Point:
        """A point on the surface (inward of the first exterior edge)."""
        xy = self._ext.coords._xy
        if xy.shape[0] < 3:
            return Point(0.0, 0.0)
        mid = 0.5 * (xy[0] + xy[1])
        t = xy[1] - xy[0]
        nrm = np.array([-t[1], t[0]], float)
        ln = float(np.linalg.norm(nrm))
        if ln < 1e-18:
            return self.centroid
        nrm = nrm / ln
        if _signed_area(xy) < 0.0:
            nrm = -nrm
        for scale in (1e-6, 1e-4, 1e-3, 1e-2, 1e-1):
            cand = Point(float(mid[0] + scale * nrm[0]), float(mid[1] + scale * nrm[1]))
            if self.contains(cand):
                return cand
        return self.centroid

    @property
    def wkb(self) -> bytes:
        from implicit_geometry_kernel.geos.wkb import dumps

        return dumps(self)

    def contains(self, other) -> bool:
        if isinstance(other, Point):
            return _point_in_ring(other.x, other.y, self._ext.coords._xy)
        if isinstance(other, LineString):
            return all(_point_in_ring(x, y, self._ext.coords._xy) for x, y in other.coords)
        return False


class MultiPolygon:
    geom_type = "MultiPolygon"

    def __init__(self, geoms=None):
        parts = []
        for g in geoms or ():
            if isinstance(g, MultiPolygon):
                parts.extend(g.geoms)
            else:
                parts.append(g if isinstance(g, Polygon) else Polygon(g))
        self._geoms = tuple(parts)

    @property
    def geoms(self):
        return self._geoms

    @property
    def is_empty(self) -> bool:
        return not self._geoms or all(g.is_empty for g in self._geoms)

    @property
    def area(self) -> float:
        return float(sum(g.area for g in self._geoms))

    @property
    def bounds(self):
        if not self._geoms:
            return (0.0, 0.0, 0.0, 0.0)
        boxes = [g.bounds for g in self._geoms if not g.is_empty]
        if not boxes:
            return (0.0, 0.0, 0.0, 0.0)
        return (
            min(b[0] for b in boxes),
            min(b[1] for b in boxes),
            max(b[2] for b in boxes),
            max(b[3] for b in boxes),
        )

    @property
    def is_valid(self) -> bool:
        return all(g.is_valid for g in self._geoms)

    @property
    def wkb(self) -> bytes:
        from implicit_geometry_kernel.geos.wkb import dumps

        return dumps(self)


class GeometryCollection:
    geom_type = "GeometryCollection"

    def __init__(self, geoms=None):
        self._geoms = tuple(geoms or ())

    @property
    def geoms(self):
        return self._geoms

    @property
    def is_empty(self) -> bool:
        return not self._geoms or all(getattr(g, "is_empty", False) for g in self._geoms)

    @property
    def area(self) -> float:
        return float(sum(float(getattr(g, "area", 0.0) or 0.0) for g in self._geoms))

    @property
    def wkb(self) -> bytes:
        from implicit_geometry_kernel.geos.wkb import dumps

        return dumps(self)


def _point_in_ring(x: float, y: float, ring: np.ndarray) -> bool:
    if ring.shape[0] < 3:
        return False
    inside = False
    x0, y0 = ring[-1]
    for x1, y1 in ring:
        if (y1 > y) != (y0 > y) and x < (x0 - x1) * (y - y1) / ((y0 - y1) + 1e-30) + x1:
            inside = not inside
        x0, y0 = x1, y1
    return inside


__all__ = [
    "CoordinateSequence",
    "GeometryCollection",
    "LineString",
    "LinearRing",
    "MultiPolygon",
    "Point",
    "Polygon",
]
