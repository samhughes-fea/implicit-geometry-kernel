"""windIO arc positions placed on the blade's own sections.

windIO gives layer edges and web attachments as ``nd_arc``: the normalised arc length
along the outer mould line, 0 at the trailing edge, round the suction side to the leading
edge and back along the pressure side to 1. The helpers in :mod:`windio_nd` map that to a
chord fraction with a fixed formula, ``x/c = 1 - 2 min(s, 1 - s)``, which treats arc
length as chord position and puts the leading edge at ``s = 0.5``; on a thick section the
two differ by several percent of chord.

Here each arc is measured on the section the blade actually has at that span (its
interpolated airfoil, unit chord), so a point at ``nd_arc = s`` is the OML point at that
arc length and its chord fraction is that point's ``x``. A blunt trailing edge's closing
segment is not part of the arc. Webs become two-foot webs (the line through the OML
points at their two arcs) and each side's spar cap keeps its own chord window, which
``BladeGeometry`` takes as ``web_positions`` pairs and ``cap_windows`` per side.

The blade passed in only has to carry the outer shape (stations and span tables); its
layout is not used.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np

from implicit_geometry_kernel.blade.windio_nd import (
    _anchors_by_name,
    _grid_to_z,
    _web_start_end,
    _z_of_grid,
    resolve_nd_arc_table,
    web_span_extent,
)

ArcMap = Callable[[float, float], tuple[float, str]]


def section_arc_loop(x: np.ndarray, upper: np.ndarray, lower: np.ndarray
                     ) -> tuple[np.ndarray, np.ndarray, int]:
    """The unit-chord OML as an open polyline TE (suction) -> LE -> TE (pressure).

    Returns ``(points, arc, i_le)``: the points, their arc fractions (0 at the suction
    trailing edge, 1 at the pressure one) and the leading-edge index.
    """
    x = np.asarray(x, float)
    su = np.column_stack([x[::-1], np.asarray(upper, float)[::-1]])
    pr = np.column_stack([x[1:], np.asarray(lower, float)[1:]])
    pts = np.vstack([su, pr])
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    return pts, cum / cum[-1], len(x) - 1


def arc_point(x: np.ndarray, upper: np.ndarray, lower: np.ndarray, s: float
              ) -> tuple[float, float, str]:
    """``(x, y, side)`` of the unit-chord OML point at arc fraction ``s``."""
    pts, arc, i_le = section_arc_loop(x, upper, lower)
    u = float(np.clip(s, 0.0, 1.0))
    px = float(np.interp(u, arc, pts[:, 0]))
    py = float(np.interp(u, arc, pts[:, 1]))
    return px, py, "suction" if u <= arc[i_le] else "pressure"


def section_arc_map(geo) -> ArcMap:
    """``f(z, s) -> (chord fraction, side)`` of arc fraction ``s`` on ``geo``'s section."""
    cache: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    def at(z: float, s: float) -> tuple[float, str]:
        key = round(float(z), 9)
        if key not in cache:
            *_, (x, upper, lower) = geo.interpolate(key)
            cache[key] = (np.asarray(x, float), np.asarray(upper, float),
                          np.asarray(lower, float))
        px, _py, side = arc_point(*cache[key], s)
        return px, side

    return at


def _clip_z(z: float, geo) -> float:
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    return float(np.clip(z, z0, z1))


def _grid_of_z(z: float, z_table: np.ndarray | None) -> float:
    """The file's span grid value at z (the inverse of ``z_of_s``)."""
    if z_table is None:
        return float(z)
    return float(np.interp(float(z), z_table[:, 1], z_table[:, 0]))


WebSeries = tuple[str, np.ndarray, np.ndarray, "tuple[float, float] | None"]


def web_feet_series_from_windio_structure(structure: Mapping[str, Any], geo,
                                          z_of_s: np.ndarray | Sequence | None = None,
                                          ) -> list[WebSeries] | None:
    """Per web, in file order: ``(name, z, feet, window)``.

    ``feet`` has shape ``(n, 2)``: the ``(suction, pressure)`` chord fractions of the
    web's ``start_nd_arc`` and ``end_nd_arc`` placed on ``geo``'s section at each grid
    point inside the web's span window ``window`` (z, or ``None`` for the whole blade;
    :func:`windio_nd.web_span_extent`). ``None`` when a web has no arc tables.
    """
    webs = [w for w in (structure.get("webs") or ()) if isinstance(w, Mapping)]
    if not webs:
        return None
    anchors = _anchors_by_name(structure)
    arc = section_arc_map(geo)
    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    series: list[WebSeries] = []
    for j, web in enumerate(webs):
        start, end = _web_start_end(web, anchors)
        if start is None or end is None:
            return None
        extent = web_span_extent(web, anchors)
        grid = np.unique(np.concatenate([start[0], end[0]]))
        if extent is not None:  # table points outside the span window are placeholders
            inside = grid[(grid >= extent[0] - 1e-12) & (grid <= extent[1] + 1e-12)]
            grid = inside if inside.size else grid
        zs = _grid_to_z(grid, z_table)
        feet = []
        for g, z in zip(grid, zs):
            a = arc(_clip_z(z, geo), float(np.interp(g, *start)))
            b = arc(_clip_z(z, geo), float(np.interp(g, *end)))
            su, pr = (a, b) if a[1] == "suction" else (b, a)
            feet.append((su[0], pr[0]))
        window = None if extent is None else (_z_of_grid(extent[0], z_table),
                                              _z_of_grid(extent[1], z_table))
        series.append((str(web.get("name", f"web_{j}")), zs, np.asarray(feet, float), window))
    return series


def web_feet_fn_from_windio_structure(structure: Mapping[str, Any], geo,
                                      z_of_s: np.ndarray | Sequence | None = None,
                                      ) -> Callable[[float], tuple[tuple[float, float], ...]] | None:
    """``f(z) ->`` the ``(suction, pressure)`` chord fractions of each web live at z,
    LE to TE (``BladeGeometry(web_positions=...)``).

    The file's arcs are interpolated to z (linear on its span grid, as windIO defines
    them) and placed on the section at z, so the feet are the outer-surface points at
    those arcs at every z, not only at the grid points: a chord fraction interpolated
    between grid points is off by millimetres where the section changes between them."""
    webs = [w for w in (structure.get("webs") or ()) if isinstance(w, Mapping)]
    if not webs:
        return None
    anchors = _anchors_by_name(structure)
    arc = section_arc_map(geo)
    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    rows = []
    for web in webs:
        start, end = _web_start_end(web, anchors)
        if start is None or end is None:
            return None
        rows.append((start, end, web_span_extent(web, anchors)))

    def at(z: float) -> tuple[tuple[float, float], ...]:
        g = _grid_of_z(float(z), z_table)
        zc = _clip_z(float(z), geo)
        out = []
        for start, end, extent in rows:
            if extent is not None and not extent[0] - 1e-12 <= g <= extent[1] + 1e-12:
                continue
            a = arc(zc, float(np.interp(g, *start)))
            b = arc(zc, float(np.interp(g, *end)))
            su, pr = (a, b) if a[1] == "suction" else (b, a)
            if 0.0 < su[0] < 1.0 and 0.0 < pr[0] < 1.0:
                out.append((su[0], pr[0]))
        out.sort(key=lambda p: p[0] + p[1])
        return tuple(out)

    return at


def _is_cap_layer(layer: Mapping[str, Any]) -> bool:
    name = str(layer.get("name", "")).lower()
    return layer.get("web") in (None, "") and any(tag in name for tag in ("spar", "cap"))


CapSeries = tuple[str, np.ndarray, np.ndarray]


def cap_side_series_from_windio_structure(structure: Mapping[str, Any], geo,
                                          z_of_s: np.ndarray | Sequence | None = None,
                                          ) -> list[CapSeries]:
    """Per spar-cap layer (by name, as :func:`windio.cap_window_table_from_windio_structure`):
    ``(side, z, rows)`` with ``rows`` the ``(lo, hi)`` chord fractions of the layer's arcs
    placed on ``geo``'s section at each grid point. A band that wraps the leading edge
    is not a spar cap and is left out."""
    anchors = _anchors_by_name(structure)
    arc = section_arc_map(geo)
    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    out: list[CapSeries] = []
    for layer in structure.get("layers") or ():
        if not isinstance(layer, Mapping) or not _is_cap_layer(layer):
            continue
        start = resolve_nd_arc_table(layer.get("start_nd_arc"), anchors)
        end = resolve_nd_arc_table(layer.get("end_nd_arc"), anchors)
        if start is None or end is None:
            continue
        grid = np.unique(np.concatenate([start[0], end[0]]))
        zs = _grid_to_z(grid, z_table)
        rows: list[tuple[float, float]] = []
        sides: set[str] = set()
        for g, z in zip(grid, zs):
            zc = _clip_z(z, geo)
            p = arc(zc, float(np.interp(g, *start)))
            q = arc(zc, float(np.interp(g, *end)))
            rows.append((min(p[0], q[0]), max(p[0], q[0])))
            sides.update((p[1], q[1]))
        if len(sides) == 1:
            out.append((sides.pop(), zs, np.asarray(rows, float)))
    return out


def _side_windows_at(series: Sequence[CapSeries], z: float) -> dict[str, tuple[float, float]]:
    out: dict[str, tuple[float, float]] = {}
    for side, zs, rows in series:
        if not zs[0] - 1e-12 <= z <= zs[-1] + 1e-12:
            continue
        lo = float(np.interp(z, zs, rows[:, 0]))
        hi = float(np.interp(z, zs, rows[:, 1]))
        if hi - lo <= 1e-6:
            continue
        if side in out:
            lo, hi = min(lo, out[side][0]), max(hi, out[side][1])
        out[side] = (lo, hi)
    return out


def cap_side_windows_fn_from_windio_structure(
    structure: Mapping[str, Any], geo, z_of_s: np.ndarray | Sequence | None = None,
) -> Callable[[float], dict[str, tuple[float, float]] | None] | None:
    """``f(z) -> {"suction": (lo, hi), "pressure": (lo, hi)}`` (``BladeGeometry
    (cap_windows=...)``): each side's cap window, the union of that side's cap layers
    live at z, their arcs interpolated to z and placed on the section at z (as
    :func:`web_feet_fn_from_windio_structure`); ``None`` where a side has none, and
    ``None`` overall when the file has no caps on both sides."""
    series = cap_side_series_from_windio_structure(structure, geo, z_of_s)
    if not {s for s, *_ in series} >= {"suction", "pressure"}:
        return None
    anchors = _anchors_by_name(structure)
    arc = section_arc_map(geo)
    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    rows = []
    for layer in structure.get("layers") or ():
        if not isinstance(layer, Mapping) or not _is_cap_layer(layer):
            continue
        start = resolve_nd_arc_table(layer.get("start_nd_arc"), anchors)
        end = resolve_nd_arc_table(layer.get("end_nd_arc"), anchors)
        if start is None or end is None:
            continue
        grid = np.unique(np.concatenate([start[0], end[0]]))
        sides = {arc(_clip_z(zg, geo), float(np.interp(g, *t)))[1]
                 for g, zg in zip(grid, _grid_to_z(grid, z_table)) for t in (start, end)}
        if len(sides) == 1:  # as the series: a band round the leading edge is not a cap
            rows.append((start, end, float(grid[0]), float(grid[-1])))

    def at(z: float) -> dict[str, tuple[float, float]] | None:
        g = _grid_of_z(float(z), z_table)
        zc = _clip_z(float(z), geo)
        out: dict[str, tuple[float, float]] = {}
        for start, end, g0, g1 in rows:
            if not g0 - 1e-12 <= g <= g1 + 1e-12:
                continue
            p = arc(zc, float(np.interp(g, *start)))
            q = arc(zc, float(np.interp(g, *end)))
            if p[1] != q[1]:
                continue
            lo, hi = min(p[0], q[0]), max(p[0], q[0])
            if hi - lo <= 1e-6:
                continue
            if p[1] in out:
                lo, hi = min(lo, out[p[1]][0]), max(hi, out[p[1]][1])
            out[p[1]] = (lo, hi)
        return out if set(out) == {"suction", "pressure"} else None

    return at


def cap_union_table_from_windio_structure(
    structure: Mapping[str, Any], geo, z_of_s: np.ndarray | Sequence | None = None,
) -> tuple[np.ndarray, np.ndarray] | None:
    """``(z, values)``, values ``(n, 2)``: the union of both sides' cap windows at each
    cap grid point where both sides have a cap (the wall-region graph's cap anchors)."""
    series = cap_side_series_from_windio_structure(structure, geo, z_of_s)
    if not series:
        return None
    zs = np.unique(np.concatenate([z for _s, z, _r in series]))
    keep, vals = [], []
    for z in zs:
        w = _side_windows_at(series, float(z))
        if set(w) == {"suction", "pressure"}:
            keep.append(float(z))
            vals.append((min(w["suction"][0], w["pressure"][0]),
                         max(w["suction"][1], w["pressure"][1])))
    if not keep:
        return None
    if len(keep) == 1:
        keep.append(keep[0] + 1e-9)
        vals.append(vals[0])
    return np.asarray(keep), np.asarray(vals, float)


__all__ = [
    "arc_point",
    "cap_side_series_from_windio_structure",
    "cap_side_windows_fn_from_windio_structure",
    "cap_union_table_from_windio_structure",
    "section_arc_loop",
    "section_arc_map",
    "web_feet_fn_from_windio_structure",
    "web_feet_series_from_windio_structure",
]
