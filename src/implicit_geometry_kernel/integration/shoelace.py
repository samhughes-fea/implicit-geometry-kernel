# src/implicit_geometry_kernel/shoelace_integration/module.py
"""
Shoelace (polygon) formulae for section properties.

Computes area, centroid, and second moments of area from polygon geometry
using the shoelace formula (Green's theorem). S frame (pivot-referenced)
and S' frame (centroid-referenced). Supports Polygon with holes and MultiPolygon.
"""

from typing import Tuple
import logging
import numpy as np

logger = logging.getLogger(__name__)


def ring_sums_S(coords_S: np.ndarray) -> Tuple[float, ...]:
    """
    Compute raw polygon sums for S frame coordinates (shoelace formula).

    Computes the raw sums required for area, centroid, and second moments
    computation using the shoelace formula (Green's theorem).

    Parameters
    ----------
    coords_S : np.ndarray
        Vertex coordinates (x, y) [m] in S frame, shape (n, 2).
        First point may or may not be repeated at end.

    Returns
    -------
    tuple of float
        (Σcross, Σ(x+x')cross, Σ(y+y')cross, Σ(y²+yy'+y'²)cross,
         Σ(x²+xx'+x'²)cross, Σ(2xy+xy'+x'y)cross).
    """
    x = coords_S[:, 0]
    y = coords_S[:, 1]
    if not np.allclose(coords_S[0], coords_S[-1]):
        x = np.append(x, x[0])
        y = np.append(y, y[0])
    cross = x[:-1] * y[1:] - x[1:] * y[:-1]
    ixy_sum = np.sum((x[:-1] * y[1:] + 2.0 * x[:-1] * y[:-1] + 2.0 * x[1:] * y[1:] + x[1:] * y[:-1]) * cross)
    result = (
        float(np.sum(cross)),
        float(np.sum((x[:-1] + x[1:]) * cross)),
        float(np.sum((y[:-1] + y[1:]) * cross)),
        float(np.sum((y[:-1] ** 2 + y[:-1] * y[1:] + y[1:] ** 2) * cross)),
        float(np.sum((x[:-1] ** 2 + x[:-1] * x[1:] + x[1:] ** 2) * cross)),
        float(ixy_sum),
    )
    if any(isinstance(v, (tuple, list)) for v in result):
        logger.error(
            "[ring_sums_S] ERROR: Returned tuple/list values! result=%s, types=%s",
            result, tuple(type(v).__name__ for v in result)
        )
        raise TypeError(f"ring_sums_S returned non-scalar values: {result}")
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug(
            "[ring_sums_S] Returning: types=(%s, %s, %s, %s, %s, %s), values=(%.6e, %.6e, %.6e, %.6e, %.6e, %.6e)",
            type(result[0]).__name__, type(result[1]).__name__, type(result[2]).__name__,
            type(result[3]).__name__, type(result[4]).__name__, type(result[5]).__name__,
            result[0], result[1], result[2], result[3], result[4], result[5]
        )
    return result


def validate_and_orient_polygon(poly, allow_empty: bool = False):
    """Validate and CCW-orient. Repair via buffer(0). Duck-typed (no Shapely import)."""
    if allow_empty and (poly.is_empty or poly.area <= 0.0):
        return poly
    if not allow_empty and (poly.is_empty or poly.area <= 0.0):
        raise ValueError("Invalid polygon: empty or zero area")
    if not poly.is_valid:
        poly = poly.buffer(0)
        if getattr(poly, "geom_type", "") == "MultiPolygon":
            raise ValueError(
                f"Polygon repair via buffer(0) produced MultiPolygon with {len(poly.geoms)} disconnected pieces."
            )
        if allow_empty and poly.is_empty:
            return poly
    if not poly.is_valid or (not allow_empty and poly.is_empty):
        raise ValueError("Invalid polygon after repair")
    return poly


def section_properties_S(
    poly_S,
) -> Tuple[float, Tuple[float, float], float, float, float]:
    """
    Compute section properties in S frame (pivot-referenced).

    Area, centroid, and second moments of area for a polygon or MultiPolygon
    in the Station (S) frame (origin at pivot, 0.333 chord). Uses shoelace formula.
    """
    if getattr(poly_S, "geom_type", "") == "MultiPolygon":
        if poly_S.is_empty:
            raise ValueError("MultiPolygon is empty")

        polygon_props = []
        for i, geom in enumerate(poly_S.geoms):
            if getattr(geom, "geom_type", "") != "Polygon":
                logger.warning("[section_properties_S] MultiPolygon contains non-Polygon at index %d, skipping", i)
                continue
            if geom.is_empty or geom.area <= 0:
                logger.warning("[section_properties_S] MultiPolygon contains empty polygon at index %d, skipping", i)
                continue

            A_i, (Cx_i, Cy_i), Ixx_i, Iyy_i, Ixy_i = section_properties_S(geom)
            polygon_props.append({
                'A': A_i, 'Cx': Cx_i, 'Cy': Cy_i,
                'Ixx': Ixx_i, 'Iyy': Iyy_i, 'Ixy': Ixy_i,
            })

        if len(polygon_props) == 0:
            raise ValueError("MultiPolygon contains no valid polygons")

        A_total = sum(p['A'] for p in polygon_props)
        if A_total <= 0:
            raise ValueError(f"Total area of MultiPolygon is non-positive: {A_total}")

        Cx_total = sum(p['A'] * p['Cx'] for p in polygon_props) / A_total
        Cy_total = sum(p['A'] * p['Cy'] for p in polygon_props) / A_total

        Ixx_total = sum(p['Ixx'] for p in polygon_props) - sum(p['A'] * (p['Cy'])**2 for p in polygon_props) + A_total * (Cy_total)**2
        Iyy_total = sum(p['Iyy'] for p in polygon_props) - sum(p['A'] * (p['Cx'])**2 for p in polygon_props) + A_total * (Cx_total)**2
        Ixy_total = sum(p['Ixy'] for p in polygon_props) - sum(p['A'] * p['Cx'] * p['Cy'] for p in polygon_props) + A_total * Cx_total * Cy_total

        return (float(A_total), (float(Cx_total), float(Cy_total)), float(Ixx_total), float(Iyy_total), float(Ixy_total))

    # Single polygon
    poly_S = validate_and_orient_polygon(poly_S, allow_empty=False)
    ring_result = ring_sums_S(np.asarray(poly_S.exterior.coords, float))
    s_cross, s_Cx, s_Cy, s_Ixx, s_Iyy, s_Ixy = ring_result

    if any(isinstance(v, (tuple, list, np.ndarray)) for v in [s_cross, s_Cx, s_Cy, s_Ixx, s_Iyy, s_Ixy]):
        raise TypeError(f"ring_sums_S (exterior) returned non-scalar: s_cross={type(s_cross)}, s_Iyy={type(s_Iyy)}, s_Ixy={type(s_Ixy)}")

    s_cross = float(s_cross)
    s_Cx = float(s_Cx)
    s_Cy = float(s_Cy)
    s_Ixx = float(s_Ixx)
    s_Iyy = float(s_Iyy)
    s_Ixy = float(s_Ixy)

    for i, hole in enumerate(poly_S.interiors):
        hole_coords = np.asarray(hole.coords, float)
        hole_result = ring_sums_S(hole_coords)
        c_cross, c_Cx, c_Cy, c_Ixx, c_Iyy, c_Ixy = hole_result

        if any(isinstance(v, (tuple, list, np.ndarray)) for v in [c_cross, c_Cx, c_Cy, c_Ixx, c_Iyy, c_Ixy]):
            raise TypeError(f"ring_sums_S (hole {i+1}) returned non-scalar")

        c_cross = float(c_cross)
        c_Cx = float(c_Cx)
        c_Cy = float(c_Cy)
        c_Ixx = float(c_Ixx)
        c_Iyy = float(c_Iyy)
        c_Ixy = float(c_Ixy)

        if c_cross < 0:
            c_cross = -c_cross
            c_Cx = -c_Cx
            c_Cy = -c_Cy
            c_Ixx = -c_Ixx
            c_Iyy = -c_Iyy
            c_Ixy = -c_Ixy

        s_cross -= c_cross
        s_Cx -= c_Cx
        s_Cy -= c_Cy
        s_Ixx -= c_Ixx
        s_Iyy -= c_Iyy
        s_Ixy -= c_Ixy

    A_S = 0.5 * s_cross
    if A_S <= 0:
        A_S = abs(A_S)

    if A_S <= 0.0:
        raise ValueError("Degenerate polygon: zero area in shoelace section_properties_S")

    Cx_S = s_Cx / (6.0 * A_S)
    Cy_S = s_Cy / (6.0 * A_S)
    Ixx_S = (1.0 / 12.0) * s_Ixx
    Iyy_S = (1.0 / 12.0) * s_Iyy
    Ixy_S = (1.0 / 24.0) * s_Ixy

    return (float(A_S), (float(Cx_S), float(Cy_S)), float(Ixx_S), float(Iyy_S), float(Ixy_S))


def section_properties_Sprime(
    poly_S,
) -> Tuple[float, Tuple[float, float], float, float, float]:
    """Section properties in S' (centroid frame) via the parallel-axis theorem."""
    A_S, (Cx_S, Cy_S), Ixx_S, Iyy_S, Ixy_S = section_properties_S(poly_S)
    Ixx_p = Ixx_S - A_S * Cy_S * Cy_S
    Iyy_p = Iyy_S - A_S * Cx_S * Cx_S
    Ixy_p = Ixy_S - A_S * Cx_S * Cy_S
    return (float(A_S), (0.0, 0.0), float(Ixx_p), float(Iyy_p), float(Ixy_p))
