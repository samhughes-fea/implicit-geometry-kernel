"""Shared inward-offset locus (G-MID). G-POLY may import this leaf."""

from __future__ import annotations

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

from implicit_geometry_kernel.contract.errors import InfeasibleSectionGeometry
from implicit_geometry_kernel.integration.shoelace import validate_and_orient_polygon

OFFSET_QUAD_SEGS = 16
OFFSET_JOIN_STYLE = 2  # mitre
OFFSET_MITRE_LIMIT = 5.0


def offset_polygon_inward(poly, depth: float, *, quad_segs: int | None = None):
    """Inward normal offset of ``poly`` by ``depth``, the one offset operator.

    Prefer this over ``poly.buffer(-d, ...)``: the join settings are the
    package's parallelism contract (bug_log SG-001 / SG-003).
    """
    d = float(depth)
    if not np.isfinite(d) or d < 0.0:
        raise ValueError(f"offset_polygon_inward requires a finite depth >= 0; got {depth!r}.")
    return poly.buffer(
        -d,
        quad_segs=OFFSET_QUAD_SEGS if quad_segs is None else int(quad_segs),
        join_style=OFFSET_JOIN_STYLE,
        mitre_limit=OFFSET_MITRE_LIMIT,
    )


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


__all__ = [
    "OFFSET_JOIN_STYLE",
    "OFFSET_MITRE_LIMIT",
    "OFFSET_QUAD_SEGS",
    "airfoil_offset_ring",
    "airfoil_polygon",
    "locus_at_depth",
    "offset_polygon_inward",
]
