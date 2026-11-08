"""PolygonOps protocol: constructive Booleans and offsets behind one seam.

Shapely/GEOS remains the overlay engine until ``KernelPolygonOps`` offset and
Boolean methods match the oracle tests. Affine helpers here are numpy. Do not
copy Shapely's Python wrappers; the algorithms are JTS/GEOS buffer (mitre) and
overlay. This module must not import Shapely (gateway is ``polygons.py``).
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable

import numpy as np

OFFSET_QUAD_SEGS = 16
OFFSET_JOIN_STYLE = 2  # mitre
OFFSET_MITRE_LIMIT = 5.0


@runtime_checkable
class PolygonOps(Protocol):
    def offset_inward(self, poly: Any, depth: float, *, quad_segs: int | None = None) -> Any: ...
    def intersection(self, a: Any, b: Any) -> Any: ...
    def difference(self, a: Any, b: Any) -> Any: ...
    def union(self, geoms) -> Any: ...
    def make_valid(self, geom: Any) -> Any: ...
    def rotate(self, geom: Any, angle_rad: float, *, cx: float = 0.0, cy: float = 0.0) -> Any: ...
    def envelope(self, minx: float, miny: float, maxx: float, maxy: float) -> Any: ...
    def orient(self, poly: Any, *, sign: float = 1.0) -> Any: ...
    def translate(self, geom: Any, xoff: float, yoff: float) -> Any: ...


def _signed_area(xy: np.ndarray) -> float:
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def oriented_rings(geom) -> list[np.ndarray]:
    """Exterior then holes as numpy rings (CCW exterior). Plotting uses this."""
    geoms = getattr(geom, "geoms", None)
    if geoms is not None:
        out: list[np.ndarray] = []
        for part in geoms:
            out.extend(oriented_rings(part))
        return out
    ext = np.asarray(geom.exterior.coords, float)
    if _signed_area(ext) < 0.0:
        ext = ext[::-1].copy()
    rings = [ext]
    rings.extend(np.asarray(r.coords, float) for r in geom.interiors)
    return rings


__all__ = [
    "OFFSET_JOIN_STYLE",
    "OFFSET_MITRE_LIMIT",
    "OFFSET_QUAD_SEGS",
    "PolygonOps",
    "oriented_rings",
]
