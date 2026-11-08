"""Little-endian ISO WKB for 2-D kernel geometries. No Shapely."""

from __future__ import annotations

import struct

import numpy as np

NDR = 1
WKB_POINT = 1
WKB_LINESTRING = 2
WKB_POLYGON = 3
WKB_MULTIPOINT = 4
WKB_MULTILINESTRING = 5
WKB_MULTIPOLYGON = 6
WKB_GEOMETRYCOLLECTION = 7
_SRID = 0x20000000
_Z = 0x80000000
_M = 0x40000000


def _close(xy: np.ndarray) -> np.ndarray:
    xy = np.asarray(xy, float)
    if xy.size == 0:
        return xy.reshape(0, 2)
    if xy.shape[0] >= 2 and np.array_equal(xy[0], xy[-1]):
        return xy
    return np.vstack([xy, xy[0]])


def _pack_xy(xy: np.ndarray) -> bytes:
    arr = np.ascontiguousarray(xy, dtype="<f8")
    return struct.pack("<I", int(arr.shape[0])) + arr.tobytes()


def dumps(geom) -> bytes:
    kind = getattr(geom, "geom_type", None)
    if kind == "Point":
        return bytes([NDR]) + struct.pack("<Idd", WKB_POINT, float(geom.x), float(geom.y))
    if kind == "LineString":
        xy = np.asarray(geom.coords, float)
        return bytes([NDR]) + struct.pack("<I", WKB_LINESTRING) + _pack_xy(xy)
    if kind == "Polygon":
        if geom.is_empty:
            return bytes([NDR]) + struct.pack("<II", WKB_POLYGON, 0)
        rings = [_close(np.asarray(geom.exterior.coords, float))]
        rings.extend(_close(np.asarray(r.coords, float)) for r in geom.interiors)
        body = struct.pack("<I", len(rings)) + b"".join(_pack_xy(r) for r in rings)
        return bytes([NDR]) + struct.pack("<I", WKB_POLYGON) + body
    if kind == "MultiPolygon":
        parts = [dumps(g) for g in geom.geoms]
        return bytes([NDR]) + struct.pack("<II", WKB_MULTIPOLYGON, len(parts)) + b"".join(parts)
    if kind == "GeometryCollection":
        parts = [dumps(g) for g in geom.geoms]
        return bytes([NDR]) + struct.pack("<II", WKB_GEOMETRYCOLLECTION, len(parts)) + b"".join(parts)
    raise TypeError(f"Cannot encode geom_type {kind!r}")


class _Reader:
    def __init__(self, data: bytes):
        self.data = data
        self.i = 0

    def u8(self) -> int:
        v = self.data[self.i]
        self.i += 1
        return v

    def take(self, n: int) -> bytes:
        sl = self.data[self.i : self.i + n]
        self.i += n
        return sl

    def u32(self, le: bool) -> int:
        fmt = "<I" if le else ">I"
        return struct.unpack(fmt, self.take(4))[0]

    def f64(self, le: bool) -> float:
        fmt = "<d" if le else ">d"
        return struct.unpack(fmt, self.take(8))[0]

    def xy(self, n: int, le: bool, dims: int) -> np.ndarray:
        pts = np.empty((n, 2), float)
        for k in range(n):
            pts[k, 0] = self.f64(le)
            pts[k, 1] = self.f64(le)
            for _ in range(dims - 2):
                self.f64(le)
        return pts


def _read_one(r: _Reader):
    from implicit_geometry_kernel.geos.geom import (
        GeometryCollection,
        LineString,
        MultiPolygon,
        Point,
        Polygon,
    )

    le = r.u8() == NDR
    raw = r.u32(le)
    typ = raw & 0xFF
    dims = 3 if raw & _Z else 2
    if raw & _M:
        dims += 1
    if raw & _SRID:
        r.u32(le)
    if typ == WKB_POINT:
        x, y = r.f64(le), r.f64(le)
        for _ in range(dims - 2):
            r.f64(le)
        return Point(x, y)
    if typ == WKB_LINESTRING:
        n = r.u32(le)
        return LineString(r.xy(n, le, dims))
    if typ == WKB_POLYGON:
        n_rings = r.u32(le)
        if n_rings == 0:
            return Polygon()
        rings = [r.xy(r.u32(le), le, dims) for _ in range(n_rings)]
        return Polygon(rings[0], rings[1:])
    if typ == WKB_MULTIPOLYGON:
        n = r.u32(le)
        return MultiPolygon([_read_one(r) for _ in range(n)])
    if typ == WKB_GEOMETRYCOLLECTION:
        n = r.u32(le)
        return GeometryCollection([_read_one(r) for _ in range(n)])
    if typ == WKB_MULTILINESTRING:
        n = r.u32(le)
        return GeometryCollection([_read_one(r) for _ in range(n)])
    if typ == WKB_MULTIPOINT:
        n = r.u32(le)
        return GeometryCollection([_read_one(r) for _ in range(n)])
    raise ValueError(f"Unsupported WKB type {typ}")


def loads(data: bytes):
    if not data:
        from implicit_geometry_kernel.geos.geom import Polygon

        return Polygon()
    return _read_one(_Reader(bytes(data)))


__all__ = ["dumps", "loads"]
