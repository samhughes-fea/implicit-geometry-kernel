"""Read windIO ``structure`` anchors into kernel chord fractions.

windIO ``nd_arc`` is normalized airfoil arc length: 0 and 1 at the trailing
edge, ~0.5 at the leading edge. A flapwise web's two feet share

``x/c = 1 - 2 min(s, 1-s)``

(LE = 0, TE = 1).

Span tables interpolate that map at each ``z``. The median helpers stay for
G8 and for span-constant layouts. ``nd_arc`` is never rewritten as
root-anchored metres.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Callable

import numpy as np


def nd_arc_to_chord_frac(s: float) -> float:
    """Map windIO ``nd_arc`` to chord fraction from the leading edge."""
    u = float(np.clip(s, 0.0, 1.0))
    return float(1.0 - 2.0 * min(u, 1.0 - u))


def _nd_table(block: Any) -> tuple[np.ndarray, np.ndarray] | None:
    if not isinstance(block, Mapping):
        return None
    if "values" not in block:
        return None
    values = np.asarray(block["values"], float).reshape(-1)
    if values.size == 0:
        return None
    grid = block.get("grid")
    if grid is None:
        grid = np.linspace(0.0, 1.0, values.size)
    return np.asarray(grid, float).reshape(-1), values


def resolve_nd_arc_table(block: Any, anchors: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray] | None:
    """Inline ``grid``/``values``, or ``anchor.name`` + ``handle`` into ``anchors``."""
    direct = _nd_table(block)
    if direct is not None:
        return direct
    if not isinstance(block, Mapping):
        return None
    ref = block.get("anchor")
    if not isinstance(ref, Mapping):
        return None
    host = anchors.get(str(ref.get("name", "")))
    if not isinstance(host, Mapping):
        return None
    handle = str(ref.get("handle", "start_nd_arc"))
    return _nd_table(host.get(handle))


def _frac_along(start: tuple[np.ndarray, np.ndarray],
                end: tuple[np.ndarray, np.ndarray] | None) -> tuple[np.ndarray, np.ndarray]:
    g, v = start
    fracs = np.empty(v.size, float)
    for i, s in enumerate(v):
        xs = [nd_arc_to_chord_frac(float(s))]
        if end is not None:
            xs.append(nd_arc_to_chord_frac(float(np.interp(float(g[i]), end[0], end[1]))))
        fracs[i] = float(np.mean(xs))
    return g, fracs


def _median_of_fracs(fracs: np.ndarray) -> float | None:
    arr = np.asarray(fracs, float)
    interior = arr[(arr > 0.02) & (arr < 0.98)]
    if interior.size == 0:
        interior = arr
    if interior.size == 0:
        return None
    return round(float(np.median(interior)), 6)


def _median_chord_frac(start: tuple[np.ndarray, np.ndarray],
                       end: tuple[np.ndarray, np.ndarray] | None) -> float | None:
    _, fracs = _frac_along(start, end)
    return _median_of_fracs(fracs)


def _anchors_by_name(structure: Mapping[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for row in structure.get("anchors") or ():
        if isinstance(row, Mapping) and row.get("name"):
            out[str(row["name"])] = row
    return out


def _web_start_end(web: Mapping[str, Any], anchors: Mapping[str, Any]):
    name = str(web.get("name", ""))
    host = anchors.get(name, web)
    start = resolve_nd_arc_table(host.get("start_nd_arc") if isinstance(host, Mapping) else None, anchors)
    if start is None:
        start = resolve_nd_arc_table(web.get("start_nd_arc"), anchors)
    end = resolve_nd_arc_table(host.get("end_nd_arc") if isinstance(host, Mapping) else None, anchors)
    if end is None:
        end = resolve_nd_arc_table(web.get("end_nd_arc"), anchors)
    return start, end


def _span_extent(block: Mapping[str, Any]) -> tuple[float, float] | None:
    """Optional ``start_nd_span`` / ``end_nd_span`` as a closed interval in grid units."""
    lo = block.get("start_nd_span", block.get("span_start"))
    hi = block.get("end_nd_span", block.get("span_end"))
    if lo is None or hi is None:
        return None
    a, b = float(lo), float(hi)
    return (min(a, b), max(a, b))


def _grid_range(table: tuple[np.ndarray, np.ndarray] | None) -> tuple[float, float] | None:
    if table is None or len(table[0]) == 0:
        return None
    return float(np.min(table[0])), float(np.max(table[0]))


def _offset_grid_range(block: Any) -> tuple[float, float] | None:
    if not isinstance(block, Mapping):
        return None
    plane = block.get("plane_intersection")
    offset = plane.get("offset") if isinstance(plane, Mapping) else None
    grid = offset.get("grid") if isinstance(offset, Mapping) else None
    if grid is None or len(grid) == 0:
        return None
    g = np.asarray(grid, float)
    return float(g.min()), float(g.max())


def web_span_extent(web: Mapping[str, Any], anchors: Mapping[str, Any]) -> tuple[float, float] | None:
    """Span window of a web in windIO grid units, or ``None`` for the whole blade.

    An explicit ``start_nd_span`` / ``end_nd_span`` wins. Otherwise the web exists
    where it is defined: the grid of its ``nd_arc`` tables, narrowed to the grid of
    its ``plane_intersection`` offset when it has one. Outside the offset grid a
    precomputed ``nd_arc`` table only holds placeholders (in IEA-10-198 every web
    reads 0.217 / 0.717 at the root), so those points are not web positions.
    """
    host = anchors.get(str(web.get("name", "")), web)
    for block in (host, web):
        if isinstance(block, Mapping):
            explicit = _span_extent(block)
            if explicit is not None:
                return explicit
    start, end = _web_start_end(web, anchors)
    ranges = [r for r in (_grid_range(start), _grid_range(end)) if r is not None]
    for block in (host, web):
        off = _offset_grid_range(block)
        if off is not None:
            ranges.append(off)
            break
    if not ranges:
        return None
    lo = max(r[0] for r in ranges)
    hi = min(r[1] for r in ranges)
    if hi < lo:
        return None
    return lo, hi


def _crop_to(grid: np.ndarray, fr: np.ndarray, extent: tuple[float, float] | None) -> tuple[np.ndarray, np.ndarray]:
    """Table points inside the web's span window (placeholders outside it dropped)."""
    if extent is None:
        return grid, fr
    keep = (grid >= extent[0] - 1e-12) & (grid <= extent[1] + 1e-12)
    if int(keep.sum()) < 1:
        return grid, fr
    return grid[keep], fr[keep]


def _grid_to_z(grid: np.ndarray, z_of_s: np.ndarray | None) -> np.ndarray:
    g = np.asarray(grid, float).reshape(-1)
    if z_of_s is None:
        return g
    table = np.asarray(z_of_s, float)
    if table.ndim != 2 or table.shape[1] < 2:
        return g
    return np.interp(g, table[:, 0], table[:, 1])


def _z_of_grid(grid_s: float, z_of_s: np.ndarray | None) -> float:
    if z_of_s is None:
        return float(grid_s)
    table = np.asarray(z_of_s, float)
    return float(np.interp(float(grid_s), table[:, 0], table[:, 1]))


def web_positions_from_windio_structure(structure: Mapping[str, Any]) -> tuple[float, ...] | None:
    """LE→TE span-median chord fractions for ``structure.webs``, or ``None`` if unusable."""
    webs = [w for w in (structure.get("webs") or ()) if isinstance(w, Mapping)]
    if not webs:
        return None
    anchors = _anchors_by_name(structure)
    fracs: list[float] = []
    for web in webs:
        start, end = _web_start_end(web, anchors)
        if start is None:
            return None
        frac = _median_chord_frac(start, end)
        if frac is None:
            return None
        fracs.append(frac)
    ordered = tuple(sorted(fracs))
    if len(ordered) != len(set(np.round(ordered, 6))) or any(not 0.0 < f < 1.0 for f in ordered):
        return None
    if any(b - a <= 0.0 for a, b in zip(ordered, ordered[1:])):
        return None
    return ordered


def web_position_table_from_windio_structure(
    structure: Mapping[str, Any],
    z_of_s: np.ndarray | Sequence | None = None,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Span table ``(z, values)`` with ``values`` shape ``(n, n_webs)``, LE→TE columns.

    Grid on each nd_arc table is windIO spanwise position; ``z_of_s`` is
    ``(s, z)`` as in example 08. Webs outside an optional span window are
    still present in the table (clamped); the callable interpolator drops them.
    """
    webs = [w for w in (structure.get("webs") or ()) if isinstance(w, Mapping)]
    if not webs:
        return None
    anchors = _anchors_by_name(structure)
    series: list[tuple[np.ndarray, np.ndarray]] = []
    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    for web in webs:
        start, end = _web_start_end(web, anchors)
        if start is None:
            return None
        g, fr = _crop_to(*_frac_along(start, end), web_span_extent(web, anchors))
        series.append((_grid_to_z(g, z_table), fr))
    z_all = np.unique(np.concatenate([z for z, _ in series]))
    if z_all.size < 1:
        return None
    if z_all.size == 1:
        z_all = np.array([z_all[0], z_all[0] + 1e-9], float)
        cols = np.column_stack([np.full(2, float(np.interp(z_all[0], z, fr))) for z, fr in series])
    else:
        cols = np.column_stack([np.interp(z_all, z, fr) for z, fr in series])
    order = np.argsort(np.median(cols, axis=0))
    cols = cols[:, order]
    if np.any(~np.isfinite(cols)):
        return None
    return z_all, cols


def web_positions_fn_from_windio_structure(
    structure: Mapping[str, Any],
    z_of_s: np.ndarray | Sequence | None = None,
) -> Callable[[float], tuple[float, ...]] | None:
    """``f(z) →`` increasing chord fractions of webs whose span window covers ``z``."""
    table = web_position_table_from_windio_structure(structure, z_of_s)
    if table is None:
        return None
    z_grid, cols = table
    webs = [w for w in (structure.get("webs") or ()) if isinstance(w, Mapping)]
    anchors = _anchors_by_name(structure)
    extents: list[tuple[float, float] | None] = []
    medians = np.median(cols, axis=0)
    web_order = np.argsort(medians)
    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    ordered_webs = [webs[int(i)] for i in web_order]
    for web in ordered_webs:
        extent = web_span_extent(web, anchors)
        if extent is None:
            extents.append(None)
        else:
            extents.append((_z_of_grid(extent[0], z_table), _z_of_grid(extent[1], z_table)))

    def at_z(z: float) -> tuple[float, ...]:
        row = tuple(float(np.interp(float(z), z_grid, cols[:, j])) for j in range(cols.shape[1]))
        kept = []
        for j, frac in enumerate(row):
            window = extents[j]
            if window is not None and not (window[0] - 1e-12 <= float(z) <= window[1] + 1e-12):
                continue
            if not 0.0 < frac < 1.0:
                continue
            kept.append(frac)
        kept.sort()
        return tuple(kept)

    return at_z
