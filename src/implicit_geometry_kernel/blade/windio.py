"""Cap windows and transverse members from windIO structure."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from implicit_geometry_kernel.blade.stations import TransverseMember
from implicit_geometry_kernel.blade.windio_nd import (  # noqa: F401
    _anchors_by_name,
    _frac_along,
    _grid_to_z,
    _median_of_fracs,
    _nd_table,
    _span_extent,
    _z_of_grid,
    nd_arc_to_chord_frac,
    resolve_nd_arc_table,
    web_position_table_from_windio_structure,
    web_positions_fn_from_windio_structure,
    web_positions_from_windio_structure,
)

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


__all__ = [
    "cap_chord_windows_from_windio_structure",
    "cap_window_table_from_windio_structure",
    "nd_arc_to_chord_frac",
    "transverse_members_from_windio",
    "web_position_table_from_windio_structure",
    "web_positions_fn_from_windio_structure",
    "web_positions_from_windio_structure",
    "web_span_windows_from_windio_structure",
]
