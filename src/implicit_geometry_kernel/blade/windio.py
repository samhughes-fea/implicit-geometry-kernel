"""Read windIO ``structure`` anchors into kernel chord fractions.

windIO ``nd_arc`` is normalized airfoil arc length: 0 and 1 at the trailing
edge, ~0.5 at the leading edge. A flapwise web's two feet share

``x/c = 1 - 2 min(s, 1-s)``

(LE = 0, TE = 1). ``BladeGeometry.web_positions`` is span-constant, so the
reader takes the median of that pair over the anchor grid.

Cap layers that name ``spar`` / ``cap`` yield a chord window from the same
map; the loft still places Family D caps from the webs, not from those
windows.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

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


def _median_chord_frac(start: tuple[np.ndarray, np.ndarray],
                       end: tuple[np.ndarray, np.ndarray] | None) -> float | None:
    g, v = start
    fracs = []
    for i, s in enumerate(v):
        xs = [nd_arc_to_chord_frac(float(s))]
        if end is not None:
            xs.append(nd_arc_to_chord_frac(float(np.interp(float(g[i]), end[0], end[1]))))
        fracs.append(float(np.mean(xs)))
    arr = np.asarray(fracs, float)
    interior = arr[(arr > 0.02) & (arr < 0.98)]
    if interior.size == 0:
        interior = arr
    if interior.size == 0:
        return None
    return round(float(np.median(interior)), 6)


def _anchors_by_name(structure: Mapping[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for row in structure.get("anchors") or ():
        if isinstance(row, Mapping) and row.get("name"):
            out[str(row["name"])] = row
    return out


def web_positions_from_windio_structure(structure: Mapping[str, Any]) -> tuple[float, ...] | None:
    """LE→TE chord fractions for ``structure.webs``, or ``None`` if unusable."""
    webs = [w for w in (structure.get("webs") or ()) if isinstance(w, Mapping)]
    if not webs:
        return None
    anchors = _anchors_by_name(structure)
    fracs: list[float] = []
    for web in webs:
        name = str(web.get("name", ""))
        host = anchors.get(name, web)
        start = resolve_nd_arc_table(host.get("start_nd_arc") if isinstance(host, Mapping) else None, anchors)
        if start is None:
            start = resolve_nd_arc_table(web.get("start_nd_arc"), anchors)
        if start is None:
            return None
        end = resolve_nd_arc_table(host.get("end_nd_arc") if isinstance(host, Mapping) else None, anchors)
        if end is None:
            end = resolve_nd_arc_table(web.get("end_nd_arc"), anchors)
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


def cap_chord_windows_from_windio_structure(
    structure: Mapping[str, Any],
) -> tuple[tuple[float, float], ...]:
    """Chord windows for layers whose names look like spar caps."""
    anchors = _anchors_by_name(structure)
    windows: list[tuple[float, float]] = []
    seen: set[tuple[float, float]] = set()
    for layer in structure.get("layers") or ():
        if not isinstance(layer, Mapping):
            continue
        name = str(layer.get("name", "")).lower()
        if layer.get("web") not in (None, ""):
            continue
        if not any(tag in name for tag in ("spar", "cap")):
            continue
        start = resolve_nd_arc_table(layer.get("start_nd_arc"), anchors)
        end = resolve_nd_arc_table(layer.get("end_nd_arc"), anchors)
        if start is None:
            continue
        lo = _median_chord_frac(start, None)
        hi = _median_chord_frac(end, None) if end is not None else lo
        if lo is None or hi is None:
            continue
        window = (min(lo, hi), max(lo, hi))
        if window[1] - window[0] < 0.02 or window[1] - window[0] > 0.95:
            continue
        key = (round(window[0], 6), round(window[1], 6))
        if key in seen:
            continue
        seen.add(key)
        windows.append(window)
    return tuple(windows)
