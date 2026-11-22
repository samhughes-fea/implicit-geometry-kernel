"""Shoelace section properties and the SDF CSG kernels against closed forms."""

from __future__ import annotations

import numpy as np
import pytest

from implicit_geometry_kernel.geos.geom import MultiPolygon, Polygon
from implicit_geometry_kernel.integration.shoelace import (
    section_properties_S,
    section_properties_Sprime,
    validate_and_orient_polygon,
)
from implicit_geometry_kernel.sdf.kernels import csg


def _box(x0: float, y0: float, w: float, h: float) -> list[tuple[float, float]]:
    return [(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)]


def test_rectangle_properties_about_origin_and_centroid() -> None:
    w, h, x0, y0 = 2.0, 1.0, 1.0, 0.5
    area, (cx, cy), ixx, iyy, ixy = section_properties_S(Polygon(_box(x0, y0, w, h)))
    assert area == pytest.approx(w * h)
    assert (cx, cy) == pytest.approx((x0 + w / 2, y0 + h / 2))
    # About the origin: centroidal value plus parallel-axis terms.
    assert ixx == pytest.approx(w * h**3 / 12 + area * cy**2)
    assert iyy == pytest.approx(h * w**3 / 12 + area * cx**2)
    assert ixy == pytest.approx(area * cx * cy)
    _a, _c, ixx_c, iyy_c, ixy_c = section_properties_Sprime(Polygon(_box(x0, y0, w, h)))
    assert (ixx_c, iyy_c, ixy_c) == pytest.approx((w * h**3 / 12, h * w**3 / 12, 0.0), abs=1e-12)


def test_a_hole_and_a_clockwise_shell_give_the_tube_values() -> None:
    outer, inner = _box(-1.0, -1.0, 2.0, 2.0), _box(-0.5, -0.5, 1.0, 1.0)
    tube = Polygon(outer[::-1], [inner])
    area, (cx, cy), ixx, iyy, ixy = section_properties_S(tube)
    assert area == pytest.approx(4.0 - 1.0)
    assert (cx, cy) == pytest.approx((0.0, 0.0), abs=1e-12)
    assert ixx == pytest.approx((2.0**4 - 1.0**4) / 12)
    assert iyy == pytest.approx(ixx)
    assert ixy == pytest.approx(0.0, abs=1e-12)


def test_multipolygon_sums_its_parts() -> None:
    left, right = Polygon(_box(-2.0, 0.0, 1.0, 1.0)), Polygon(_box(1.0, 0.0, 1.0, 1.0))
    area, (cx, cy), ixx, iyy, _ixy = section_properties_S(MultiPolygon([left, right]))
    parts = [section_properties_S(p) for p in (left, right)]
    assert area == pytest.approx(2.0)
    assert (cx, cy) == pytest.approx((0.0, 0.5))
    assert ixx == pytest.approx(sum(p[2] for p in parts))
    assert iyy == pytest.approx(sum(p[3] for p in parts))


def test_empty_polygons_are_refused_unless_allowed() -> None:
    with pytest.raises(ValueError):
        validate_and_orient_polygon(Polygon())
    assert validate_and_orient_polygon(Polygon(), allow_empty=True).is_empty


def _circle(cx: float, cy: float, r: float):
    return lambda x, y: np.hypot(np.asarray(x) - cx, np.asarray(y) - cy) - r


X = np.linspace(-2.0, 2.0, 41)
Y = np.zeros_like(X)


def test_boolean_kernels_are_min_max_of_their_fields() -> None:
    a, b = _circle(-0.5, 0.0, 1.0), _circle(0.5, 0.0, 1.0)
    fa, fb = a(X, Y), b(X, Y)
    np.testing.assert_allclose(csg.union(a, b)(X, Y), np.minimum(fa, fb))
    np.testing.assert_allclose(csg.intersect(a, b)(X, Y), np.maximum(fa, fb))
    np.testing.assert_allclose(csg.subtract(a, b)(X, Y), np.maximum(fa, -fb))
    np.testing.assert_allclose(csg.union_all([a, b, 0.25])(X, Y), np.minimum(np.minimum(fa, fb), 0.25))
    np.testing.assert_allclose(csg.intersect_all([a, b])(X, Y), np.maximum(fa, fb))
    assert csg.union_all([a]) is a
    np.testing.assert_allclose(csg.intersect_all([0.5])(X, Y), 0.5)
    with pytest.raises(ValueError):
        csg.union_all([])
    with pytest.raises(ValueError):
        csg.intersect_all([])


def test_smooth_kernels_bound_the_sharp_ones_and_match_far_away() -> None:
    a, b = _circle(-0.5, 0.0, 1.0), _circle(0.5, 0.0, 1.0)
    fa, fb = a(X, Y), b(X, Y)
    k = 0.2
    su, si = csg.smooth_union(a, b, k)(X, Y), csg.smooth_intersect(a, b, k)(X, Y)
    assert np.all(su <= np.minimum(fa, fb) + 1e-12)
    assert np.all(si >= np.maximum(fa, fb) - 1e-12)
    far = np.abs(fa - fb) > k
    np.testing.assert_allclose(su[far], np.minimum(fa, fb)[far])
    np.testing.assert_allclose(
        csg.smooth_subtract(a, b, k)(X, Y), csg.smooth_intersect(a, lambda x, y: -b(x, y), k)(X, Y)
    )


def test_offset_shell_and_blend() -> None:
    a, b = _circle(0.0, 0.0, 1.0), _circle(0.0, 0.0, 0.5)
    fa, fb = a(X, Y), b(X, Y)
    np.testing.assert_allclose(csg.offset(a, 0.1)(X, Y), fa - 0.1)
    np.testing.assert_allclose(csg.shell(a, 0.2)(X, Y), np.abs(fa) - 0.1)
    np.testing.assert_allclose(csg.blend(a, b, 0.25)(X, Y), 0.75 * fa + 0.25 * fb)
