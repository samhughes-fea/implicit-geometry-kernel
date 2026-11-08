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
        g, fr = _frac_along(start, end)
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
        host = anchors.get(str(web.get("name", "")), web)
        extent = None
        if isinstance(host, Mapping):
            extent = _span_extent(host)
        if extent is None:
            extent = _span_extent(web)
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


def cap_chord_windows_from_windio_structure(
    structure: Mapping[str, Any],
) -> tuple[tuple[float, float], ...]:
    """Span-median chord windows for layers whose names look like spar caps."""
    rows = _cap_window_series(structure)
    windows: list[tuple[float, float]] = []
    seen: set[tuple[float, float]] = set()
    for _g, lo, hi in rows:
        window = (min(lo, hi), max(lo, hi))
        if window[1] - window[0] < 0.02 or window[1] - window[0] > 0.95:
            continue
        key = (round(window[0], 6), round(window[1], 6))
        if key in seen:
            continue
        seen.add(key)
        windows.append(window)
    return tuple(windows)


def _cap_window_series(structure: Mapping[str, Any]) -> list[tuple[np.ndarray, float, float]]:
    """Median (lo, hi) plus the start-table grid for each cap layer."""
    anchors = _anchors_by_name(structure)
    out: list[tuple[np.ndarray, float, float]] = []
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
        g, lo_s = _frac_along(start, None)
        lo = _median_of_fracs(lo_s)
        if end is not None:
            _, hi_s = _frac_along(end, None)
            hi = _median_of_fracs(hi_s)
        else:
            hi = lo
        if lo is None or hi is None:
            continue
        out.append((g, lo, hi))
    return out


def _cap_layer_frac_tables(structure: Mapping[str, Any]) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    anchors = _anchors_by_name(structure)
    out: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
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
        gs, lo = _frac_along(start, None)
        if end is None:
            hi = lo
            ge = gs
        else:
            ge, hi = _frac_along(end, None)
        g = np.unique(np.concatenate([gs, ge]))
        lo_i = np.interp(g, gs, lo)
        hi_i = np.interp(g, ge, hi)
        out.append((g, np.minimum(lo_i, hi_i), np.maximum(lo_i, hi_i)))
    return out


def cap_window_table_from_windio_structure(
    structure: Mapping[str, Any],
    z_of_s: np.ndarray | Sequence | None = None,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Union chord window as ``(z, values)`` with ``values`` shape ``(n, 2)`` (lo, hi)."""
    layers = _cap_layer_frac_tables(structure)
    if not layers:
        return None
    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    z_all = np.unique(np.concatenate([_grid_to_z(g, z_table) for g, _lo, _hi in layers]))
    if z_all.size == 1:
        z_all = np.array([z_all[0], z_all[0] + 1e-9], float)
    lo_cols = []
    hi_cols = []
    for g, lo, hi in layers:
        zg = _grid_to_z(g, z_table)
        lo_cols.append(np.interp(z_all, zg, lo))
        hi_cols.append(np.interp(z_all, zg, hi))
    lo = np.min(np.column_stack(lo_cols), axis=1)
    hi = np.max(np.column_stack(hi_cols), axis=1)
    wide = (hi - lo) >= 0.02
    if not np.any(wide):
        return None
    return z_all, np.column_stack([lo, hi])


def web_span_windows_from_windio_structure(
    structure: Mapping[str, Any],
    z_of_s: np.ndarray | Sequence | None = None,
) -> tuple[tuple[float, float] | None, ...]:
    """Per-web span windows in z, ordered LE→TE like :func:`web_position_table_from_windio_structure`."""
    table = web_position_table_from_windio_structure(structure, z_of_s)
    if table is None:
        return ()
    _z_grid, cols = table
    webs = [w for w in (structure.get("webs") or ()) if isinstance(w, Mapping)]
    anchors = _anchors_by_name(structure)
    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    order = np.argsort(np.median(cols, axis=0))
    out: list[tuple[float, float] | None] = []
    for i in order:
        web = webs[int(i)] if int(i) < len(webs) else {}
        host = anchors.get(str(web.get("name", "")), web)
        extent = _span_extent(host) if isinstance(host, Mapping) else None
        if extent is None:
            extent = _span_extent(web) if isinstance(web, Mapping) else None
        if extent is None:
            out.append(None)
        else:
            out.append((_z_of_grid(extent[0], z_table), _z_of_grid(extent[1], z_table)))
    return tuple(out)


def _layer_thickness_at_z(structure: Mapping[str, Any], z: float,
                          z_of_s: np.ndarray | Sequence | None,
                          predicate) -> float:
    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    total = 0.0
    found = False
    for layer in structure.get("layers") or ():
        if not isinstance(layer, Mapping):
            continue
        name = str(layer.get("name", ""))
        if layer.get("web") not in (None, ""):
            continue
        if not predicate(name):
            continue
        block = layer.get("thickness")
        table = _nd_table(block)
        if table is None:
            continue
        g, v = table
        zg = _grid_to_z(g, z_table)
        total += float(np.interp(float(z), zg, v))
        found = True
    return total if found else 0.012


def _is_skin_or_te_layer(name: str) -> bool:
    key = name.lower()
    if any(w in key for w in ("foam", "filler", "balsa", "sandwich")) or key.endswith("_core"):
        return False
    if key.startswith("dp") and "uniax" not in key and "triax" not in key:
        return False
    return "shell" in key or "skin" in key or "triax" in key or key in {"te", "trailing_edge"} or "te_reinforcement" in key


def _yaml_transverse_rows(structure: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    rows: list[Mapping[str, Any]] = []
    for key in ("ribs", "bulkheads", "transverse"):
        block = structure.get(key) or ()
        if isinstance(block, Mapping):
            block = (block,)
        for row in block:
            if isinstance(row, Mapping):
                rows.append(row)
    return rows


def transverse_members_from_windio(
    structure: Mapping[str, Any],
    z_of_s: np.ndarray | Sequence | None = None,
) -> tuple[TransverseMember, ...]:
    """Tip-closure plate plus any YAML ``ribs`` / ``bulkheads`` / ``transverse`` rows.

    windIO files typically omit a tip block. The plate sits at last span and takes
    the tip skin / TE gauge (geometry thickness only).
    """
    from implicit_geometry_kernel.blade.stations import TransverseMember

    z_table = None if z_of_s is None else np.asarray(z_of_s, float)
    if z_table is not None and z_table.size:
        arr = np.asarray(z_table, float)
        z_tip = float(arr[-1, 1]) if arr.ndim == 2 and arr.shape[1] >= 2 else float(arr.reshape(-1)[-1])
    else:
        z_tip = 1.0
    t_tip = _layer_thickness_at_z(structure, z_tip, z_table, _is_skin_or_te_layer)
    if t_tip <= 0.0:
        t_tip = 0.012
    members = [TransverseMember("tip_closure", float(z_tip), float(t_tip))]
    for row in _yaml_transverse_rows(structure):
        name = str(row.get("name", "rib")).strip() or "rib"
        z = row.get("z")
        if z is None:
            nd = row.get("nd_span", row.get("spanwise_position", row.get("grid")))
            if nd is None:
                continue
            z = _z_of_grid(float(np.asarray(nd, float).reshape(-1)[0]), z_table)
        thick = row.get("thickness", row.get("t"))
        if thick is None:
            block = row.get("thickness_table")
            table = _nd_table(block) if isinstance(block, Mapping) else None
            thick = float(table[1][0]) if table is not None else t_tip
        members.append(TransverseMember(name, float(z), float(thick)))
    # unique z, tip_closure kept even if a YAML rib shares the tip
    by_key: dict[tuple[str, float], TransverseMember] = {}
    for m in members:
        by_key[(m.name, round(float(m.z), 9))] = m
    return tuple(sorted(by_key.values(), key=lambda m: (m.z, m.name)))
