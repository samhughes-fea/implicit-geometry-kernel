"""ctypes libgeos_c wrapper: mitre inward offset and overlay on a box."""

from __future__ import annotations

from implicit_geometry_kernel.geos import lib as geos_lib
from implicit_geometry_kernel.geos.geom import Polygon
from implicit_geometry_kernel.geos.wkb import dumps, loads
from implicit_geometry_kernel.section.boundary.polygons import (
    OFFSET_JOIN_STYLE,
    OFFSET_MITRE_LIMIT,
    OFFSET_QUAD_SEGS,
    GeosPolygonOps,
)


def _box():
    return Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 0.4), (0.0, 0.4)])


def test_ctypes_mitre_inward_box_area():
    poly = _box()
    via_geos = loads(
        geos_lib.buffer(
            dumps(poly),
            -0.05,
            quad_segs=OFFSET_QUAD_SEGS,
            join_style=OFFSET_JOIN_STYLE,
            mitre_limit=OFFSET_MITRE_LIMIT,
        )
    )
    assert abs(via_geos.area - 0.9 * 0.3) < 1e-9


def test_ctypes_overlay_partition_box():
    poly = _box()
    clip = Polygon([(0.2, -0.1), (0.7, -0.1), (0.7, 0.5), (0.2, 0.5)])
    inter = loads(geos_lib.intersection(dumps(poly), dumps(clip)))
    diff = loads(geos_lib.difference(dumps(poly), dumps(clip)))
    assert abs(inter.area + diff.area - poly.area) < 1e-9


def test_geos_polygon_ops_is_the_production_engine():
    inner = GeosPolygonOps().offset_inward(_box(), 0.05)
    assert abs(inner.area - 0.9 * 0.3) < 1e-9
