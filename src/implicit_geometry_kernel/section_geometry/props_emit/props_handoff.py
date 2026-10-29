"""Peer handoffs from geometry → ``section_properties``.

Geometry emits outlines / void loops; props aggregates. These helpers pack the
common seam so demos and BSP do not invent fiber extents or warping topology.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike


def fiber_extents_from_outline(
    outline_xy: ArrayLike,
    *,
    Cx: float,
    Cy: float,
) -> dict[str, float] | None:
    """OML AABB half-extents about ``(Cx, Cy)`` for section moduli ``W_x, W_y``.

    Returns ``{"y_max": …, "x_max": …}`` or ``None`` when the outline is empty
    or degenerate.
    """
    pts = np.asarray(outline_xy, dtype=float)
    if pts.size == 0 or pts.ndim != 2 or pts.shape[1] != 2:
        return None
    y_max = float(np.max(np.abs(pts[:, 1] - float(Cy))))
    x_max = float(np.max(np.abs(pts[:, 0] - float(Cx))))
    if y_max <= 0.0 or x_max <= 0.0:
        return None
    return {"y_max": y_max, "x_max": x_max}


def fiber_extents_from_airfoil(
    airfoil: Any,
    *,
    Cx: float,
    Cy: float,
) -> dict[str, float] | None:
    """Fiber extents from an airfoil OML (``vertices``) about ``(Cx, Cy)``."""
    verts = getattr(airfoil, "vertices", None)
    if verts is None:
        return None
    return fiber_extents_from_outline(verts, Cx=float(Cx), Cy=float(Cy))


def outer_polygon_from_airfoil(airfoil: Any):
    """Shapely polygon of the airfoil OML in geometry/NACA frame (warping ``Jt``).

    Returns ``None`` when shapely is unavailable or the outline is invalid.
    """
    try:
        from shapely.geometry import Polygon
    except ImportError:
        return None
    verts = np.asarray(getattr(airfoil, "vertices", None), dtype=float)
    if verts is None or verts.ndim != 2 or verts.shape[0] < 3:
        return None
    try:
        poly = Polygon(verts)
        if not poly.is_valid or poly.area <= 0.0:
            poly = poly.buffer(0)
        return poly if poly.area > 0.0 else None
    except Exception:  # noqa: BLE001
        return None


def _cell_loop_as_bredt_mapping(
    loop: Any,
    *,
    thickness: float | None,
) -> Mapping[str, Any] | None:
    """Adapt a ``CellWallLoop`` (or mapping) for ``warping.bredt_jt_from_cell_loops``."""
    if thickness is None or float(thickness) <= 0.0:
        return None
    if isinstance(loop, Mapping):
        area = loop.get("area", loop.get("enclosed_area", loop.get("area_m2")))
        xy = loop.get(
            "xy",
            loop.get(
                "loop_xy",
                loop.get(
                    "points",
                    loop.get("polyline_twisted_pitch_axis_blade_frame"),
                ),
            ),
        )
        t = loop.get("thickness", loop.get("t", thickness))
        if area is None or xy is None:
            return None
        return {"area": float(area), "xy": xy, "thickness": float(t)}
    area = getattr(loop, "area_m2", None)
    xy = getattr(loop, "polyline_twisted_pitch_axis_blade_frame", None)
    if area is None or xy is None:
        return None
    return {"area": float(area), "xy": xy, "thickness": float(thickness)}


def topology_for_section_properties(
    *,
    airfoil: Any | None = None,
    geom_result: Any | None = None,
    outer_polygon: Any | None = None,
    cell_wall_loops: Sequence[Any] | None = None,
    bredt_thickness: float | None = None,
) -> dict[str, Any]:
    """Pack ``topology=`` for ``compute_station_section_properties``.

    Prefers an explicit ``outer_polygon``, else airfoil OML. Optional
    ``cell_wall_loops`` come from ``geom_result.cell_wall_loops`` (void-loops
    gate) when not passed explicitly. Bredt mappings are filled only when
    ``bredt_thickness`` is provided (loops alone lack wall thickness).
    """
    topo: dict[str, Any] = {}
    poly = outer_polygon
    if poly is None and airfoil is not None:
        poly = outer_polygon_from_airfoil(airfoil)
    if poly is not None:
        topo["outer_polygon"] = poly

    loops = cell_wall_loops
    if loops is None and geom_result is not None:
        loops = getattr(geom_result, "cell_wall_loops", None) or None
    if loops:
        if bredt_thickness is not None:
            mapped = [
                m
                for m in (
                    _cell_loop_as_bredt_mapping(loop, thickness=bredt_thickness)
                    for loop in loops
                )
                if m is not None
            ]
            if mapped:
                topo["cell_wall_loops"] = mapped
        else:
            topo["cell_wall_loops"] = tuple(loops)
    return topo


__all__ = [
    "fiber_extents_from_outline",
    "fiber_extents_from_airfoil",
    "outer_polygon_from_airfoil",
    "topology_for_section_properties",
]
