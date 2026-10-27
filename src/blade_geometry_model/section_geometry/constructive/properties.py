"""Shoelace section properties from constructive wall-solid polygons."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

from blade_geometry_model.shoelace_integration import (
    section_properties_S,
)


def _empty_props() -> dict[str, float]:
    nan = float("nan")
    return {
        "area": 0.0,
        "cx": nan,
        "cy": nan,
        "Ixx": nan,
        "Iyy": nan,
        "Ixy": nan,
        "r_gyr_x": nan,
        "r_gyr_y": nan,
    }


def _origin_moments(poly: Polygon | MultiPolygon) -> tuple[float, float, float, float, float, float]:
    """Return (A, Sx, Sy, Ixx0, Iyy0, Ixy0) with second moments about the origin."""
    if isinstance(poly, MultiPolygon):
        area = sx = sy = ixx0 = iyy0 = ixy0 = 0.0
        for g in poly.geoms:
            if g.is_empty or float(g.area) <= 0.0:
                continue
            try:
                a, (cx, cy), ixx, iyy, ixy = section_properties_S(g)
            except (ValueError, ZeroDivisionError):
                continue
            if a <= 0.0:
                continue
            # section_properties_S returns origin-referenced moments for a single Polygon.
            area += a
            sx += a * cx
            sy += a * cy
            ixx0 += ixx
            iyy0 += iyy
            ixy0 += ixy
        return area, sx, sy, ixx0, iyy0, ixy0

    a, (cx, cy), ixx, iyy, ixy = section_properties_S(poly)
    if a <= 0.0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    return a, a * cx, a * cy, ixx, iyy, ixy


def _props_from_poly(poly: Polygon | MultiPolygon) -> dict[str, float]:
    if poly is None or poly.is_empty or float(poly.area) <= 0.0:
        return _empty_props()
    try:
        area, sx, sy, ixx0, iyy0, ixy0 = _origin_moments(poly)
    except (ValueError, ZeroDivisionError):
        return _empty_props()
    if area <= 0.0:
        return _empty_props()
    cx = sx / area
    cy = sy / area
    ixx_c = float(ixx0) - float(area) * float(cy) ** 2
    iyy_c = float(iyy0) - float(area) * float(cx) ** 2
    ixy_c = float(ixy0) - float(area) * float(cx) * float(cy)
    return {
        "area": float(area),
        "cx": float(cx),
        "cy": float(cy),
        "Ixx": float(ixx_c),
        "Iyy": float(iyy_c),
        "Ixy": float(ixy_c),
        "r_gyr_x": float(np.sqrt(max(ixx_c, 0.0) / area)) if area > 0 else float("nan"),
        "r_gyr_y": float(np.sqrt(max(iyy_c, 0.0) / area)) if area > 0 else float("nan"),
    }


def section_properties_from_polygons(
    label_polys: Mapping[str, Any],
) -> dict[str, dict[str, float]]:
    """Compute per-label section properties matching the quadtree report schema."""
    out: dict[str, dict[str, float]] = {}
    for label, poly in label_polys.items():
        props = _props_from_poly(poly)
        props["label"] = str(label)
        out[str(label)] = props
    return out
