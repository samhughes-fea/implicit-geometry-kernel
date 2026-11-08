"""Kernel GEOS gateway: ctypes ``libgeos_c``, WKB, numpy-backed 2-D types.

Offset and overlay stay GEOS (JTS C++). This package is the swap surface for a
C++ translation. Set ``GEOS_LIBRARY_PATH`` or install system/conda GEOS.
"""

from implicit_geometry_kernel.geos import lib
from implicit_geometry_kernel.geos.geom import (
    GeometryCollection,
    LineString,
    MultiPolygon,
    Point,
    Polygon,
)
from implicit_geometry_kernel.geos.wkb import dumps, loads

__all__ = [
    "GeometryCollection",
    "LineString",
    "MultiPolygon",
    "Point",
    "Polygon",
    "dumps",
    "lib",
    "loads",
]
