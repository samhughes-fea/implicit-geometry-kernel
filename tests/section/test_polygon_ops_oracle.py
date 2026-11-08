"""GeosPolygonOps / KernelPolygonOps on the live operation subset."""

from __future__ import annotations

from implicit_geometry_kernel.geos.geom import Polygon
from implicit_geometry_kernel.section.boundary.polygons import GeosPolygonOps, KernelPolygonOps

_ATOL_AREA = 1e-12
_ATOL_HAUSDORFF = 1e-9


def _box():
    return Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 0.4), (0.0, 0.4)])


def test_kernel_rotate_matches_geos_ops_area_and_hausdorff():
    geos_ops = GeosPolygonOps()
    kernel = KernelPolygonOps()
    poly = _box()
    a = geos_ops.rotate(poly, 0.3, cx=0.2, cy=0.1)
    b = kernel.rotate(poly, 0.3, cx=0.2, cy=0.1)
    assert abs(a.area - b.area) < _ATOL_AREA
    assert a.exterior.hausdorff_distance(b.exterior) < _ATOL_HAUSDORFF


def test_kernel_offset_and_overlay_on_box():
    kernel = KernelPolygonOps()
    poly = _box()
    clip = Polygon([(0.2, -0.1), (0.7, -0.1), (0.7, 0.5), (0.2, 0.5)])
    inner = kernel.offset_inward(poly, 0.05)
    assert abs(inner.area - 0.9 * 0.3) < 1e-9
    inter = kernel.intersection(poly, clip)
    diff = kernel.difference(poly, clip)
    assert abs(inter.area + diff.area - poly.area) < 1e-9


def test_kernel_infeasible_empty_offset_is_empty():
    kernel = KernelPolygonOps()
    tiny = Polygon([(0.0, 0.0), (0.02, 0.0), (0.02, 0.02), (0.0, 0.02)])
    out = kernel.offset_inward(tiny, 0.05)
    assert out.is_empty or out.area == 0.0
