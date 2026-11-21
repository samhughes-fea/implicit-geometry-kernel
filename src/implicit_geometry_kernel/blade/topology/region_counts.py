"""Node counts per region from its length, at a target element size ``h``.

A block is the run of stations on which a region exists. The block takes one
edge count, so both ends of a bay carry the same polyline and the loft is
ruled quads. The count is the most edges any station asks for, from
:func:`_even_edges` of that station's arc length and its target edge. A
shorter station keeps that count, so its edges fall below ``h / 2``. A run
that opens or closes is the collapsed face of the loft: the count then stops
at the most edges that fan can carry, and the rest of the run shares it.
Counts are odd (an even number of edges per region).

A section without webs is one cell whose ring is the whole skin loop; the
loop then gets at least :data:`MIN_CELL_NODES` nodes (the longest region
takes the extra edges), so a rib or tip closure there is a polygon the
polygon mesher can fill rather than a quad it leaves out. A cell bounded by a
web already has at least six.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Mapping, Sequence

import numpy as np

from implicit_geometry_kernel.blade.topology.sizing import EPS_OVER_C, KAPPA_FLOOR, polyline_kappa
from implicit_geometry_kernel.blade.topology.template import ShellTopologyTemplate


MIN_CELL_NODES = 5
"""Fewest nodes round a rib or tip closure cell (polygons of at least five nodes)."""


def _length(pts: np.ndarray) -> float:
    pts = np.asarray(pts, float)
    if len(pts) < 2:
        return 0.0
    return float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))


def _even_edges(length: float, h: float, minimum: int = 3) -> int:
    edges = max(int(np.ceil(length / h - 1e-9)), minimum - 1)
    return edges + edges % 2


def _odd_nodes(length: float, h: float, minimum: int = 3) -> int:
    return _even_edges(length, h, minimum) + 1


def _held_edges(length: float, previous: int | None, h: float, minimum: int) -> tuple[int, str]:
    """``(edges, reason)``. Hold ``previous`` while ``length / previous`` is in ``[h / 2, h]``.

    Outside that band the count moves two edges toward the count ``length`` asks
    for (``step``). The first station, with no ``previous``, takes the length
    count (``length``).
    """
    edges = _even_edges(length, h, minimum)
    if previous is None:
        return edges, "length"
    lo, hi = 0.5 * h * previous, h * previous
    if lo - 1e-9 <= length <= hi + 1e-9:
        return previous, "held"
    if edges > previous:
        return min(edges, previous + 2), "step"
    if edges < previous:
        return max(edges, previous - 2), "step"
    return previous, "held"


_SITE_CORNER_MIN = 20.0
"""A 2:1 site whose inner corner would fall under this stays at the previous count."""


def _bay_chord_ratio(chords: Sequence[float] | None, index: int) -> float:
    """``|c1 - c0| / max(c0, c1)`` of the bay ending at ``index``. Zero when the chord is unknown."""
    if chords is None or index <= 0:
        return 0.0
    c0 = float(chords[index - 1])
    c1 = float(chords[index])
    den = max(c0, c1)
    if den <= 0.0:
        return 0.0
    return abs(c1 - c0) / den


def _site_corner_deg(edge: float, bay: float, chord_ratio: float) -> float:
    """Acute corner of a full-bay 2:1 trapezoid, in degrees.

    The slant runs from the coarse edge to the fine node across the whole bay,
    so the corner is ``atan(bay / (edge / 2))``. A chord change shortens that
    rise. A long bay, or a short edge, clears 20°.
    """
    if bay <= 1e-12 or edge <= 1e-12:
        return 90.0
    lift = float(bay) * max(0.0, 1.0 - float(chord_ratio))
    return math.degrees(math.atan2(lift, 0.5 * edge))


def _hold_sharp_sites(taken: dict[str, int], reasons: dict[str, str], incoming: dict[str, int],
                      smp: Mapping[str, np.ndarray], bay: float, chord_ratio: float,
                      targets: Mapping[str, tuple[float, str]]) -> None:
    """Revert one region's step when its new edge is on target and the site is under 20°.

    The revert is that region only. A region whose new edge is still outside
    ``[h / 2, h]`` keeps the two-edge step, including when the corner is under
    20°, so a long flat arc reaches its own target instead of freezing the
    station. A region's first station has no previous count and is left as the
    length asked.
    """
    for name, count in list(taken.items()):
        prev = incoming.get(name)
        if prev is None or count == prev or count <= 0:
            continue
        length = _length(smp[name]) if name in smp else 0.0
        edge = length / float(count)
        target = float(targets.get(name, (0.0, ""))[0])
        if target > 0.0 and (edge < 0.5 * target - 1e-9 or edge > target + 1e-9):
            continue
        if _site_corner_deg(edge, bay, chord_ratio) < _SITE_CORNER_MIN:
            taken[name] = prev
            reasons[name] = "taper"


@dataclass(frozen=True)
class CountDecision:
    """One region's edge count at one station, and why that count was taken."""

    z: float
    region: str
    length: float
    edge: float
    edge_from: str
    previous: int
    count: int
    reason: str


def _sagitta_edge(pts: np.ndarray, chord: float) -> float:
    """Edge the region's own polyline can carry at ``EPS_OVER_C`` of chord. Infinite when the curve is flat."""
    pts = np.asarray(pts, float)
    if len(pts) < 3 or chord <= 0.0:
        return math.inf
    kap = polyline_kappa(pts)
    k = float(np.max(kap)) if len(kap) else 0.0
    if k <= KAPPA_FLOOR:
        return math.inf
    return math.sqrt(8.0 * EPS_OVER_C * float(chord) / k)


def _fan_edges(length: float, bay: float) -> int:
    """Most even edges a collapsed arc of ``length`` can carry across ``bay``.

    The collapsed face splits each coarse segment into three quads. A segment
    that subtends much under about 28 degrees at the far station leaves a quad
    corner under 20 degrees, so the count stops there.
    """
    if bay <= 1e-6 or length <= 1e-9:
        return 10**6
    n_coarse = int(math.floor(length / (2.0 * bay * math.tan(math.radians(15.0))) + 1e-9))
    return max(1, n_coarse) * 2


def _block_counts(asked: Sequence[Mapping[str, int]], lengths: Sequence[Mapping[str, float]],
                  zs: Sequence[float], skin: set[str],
                  h_at: Sequence[float] | None = None) -> list[dict[str, int]]:
    """One edge count per region on each run of stations where that region exists.

    The count is the largest ask in the run. A station where the region is
    absent ends the run, so a region that dies and is born again takes a new
    count. A skin run that opens or closes keeps only as many edges as the
    collapsed face can carry, unless that would make any station's edges
    longer than its sizing ceiling ``h``. A web that stops is omitted, so it
    keeps the ask.
    """
    names: set[str] = set()
    for station in asked:
        names.update(station)
    out = [dict(station) for station in asked]
    for name in names:
        start: int | None = None
        for index, station in enumerate(asked):
            if name in station:
                if start is None:
                    start = index
                continue
            if start is not None:
                _stamp_block(out, lengths, zs, skin, name, start, index, h_at)
                start = None
        if start is not None:
            _stamp_block(out, lengths, zs, skin, name, start, len(asked), h_at)
    return out


def _stamp_block(out: list[dict[str, int]], lengths: Sequence[Mapping[str, float]],
                 zs: Sequence[float], skin: set[str], name: str, start: int, stop: int,
                 h_at: Sequence[float] | None = None) -> None:
    count = max(int(out[i][name]) for i in range(start, stop))
    if name in skin:
        budgets: list[int] = []
        if start > 0:
            budgets.append(_fan_edges(
                float(lengths[start][name]), abs(float(zs[start]) - float(zs[start - 1]))))
        if stop < len(out):
            budgets.append(_fan_edges(
                float(lengths[stop - 1][name]), abs(float(zs[stop - 1]) - float(zs[stop]))))
        if budgets:
            fan = min(budgets)
            if h_at is None:
                count = min(count, fan)
            else:
                # Accept the fan only when every station in the run still meets h.
                ok = all(
                    float(lengths[i][name]) / float(fan) <= float(h_at[i])
                    for i in range(start, stop)
                    if name in lengths[i] and fan > 0
                )
                if ok:
                    count = min(count, fan)
    for index in range(start, stop):
        out[index][name] = count


def _widen_single_cell(tmpl: ShellTopologyTemplate, smp: Mapping[str, np.ndarray],
                       held: dict[str, int], minimum: int) -> None:
    """Add two edges at a time to the longest loop region until the loop closes on
    :data:`MIN_CELL_NODES` nodes (a closed ring has as many nodes as edges)."""
    names = [b.name for b in tmpl.loop if b.name in smp]
    if not names:
        return
    for name in names:
        held.setdefault(name, _even_edges(0.0, 1.0, minimum))
    longest = max(names, key=lambda n: _length(smp[n]))
    while sum(held[n] for n in names) < MIN_CELL_NODES:
        held[longest] += 2


def _station_sizes(h: float | Sequence[float], n: int) -> list[float]:
    """One positive target edge per station. A scalar is that edge at every station."""
    arr = np.asarray(h, float)
    if arr.ndim == 0:
        if float(arr) <= 0.0:
            raise ValueError("element size h must be positive")
        return [float(arr)] * n
    if arr.shape != (n,):
        raise ValueError("one element size per station")
    if np.any(arr <= 0.0):
        raise ValueError("element size h must be positive")
    return [float(v) for v in arr]


def _region_target(name: str, pts: np.ndarray, ceiling: float, chord: float | None,
                   overrides: Mapping[str, float], webs: set[str]) -> tuple[float, str]:
    """``(edge, where it came from)``. Webs take the ceiling. Skin takes the tighter of the ceiling and the sagitta."""
    if name in overrides:
        return overrides[name], "region_h"
    if name in webs or chord is None:
        return ceiling, "ceiling"
    sag = _sagitta_edge(pts, chord)
    if sag < ceiling:
        return sag, "sagitta"
    return ceiling, "ceiling"


def region_counts(templates: Sequence[ShellTopologyTemplate],
                  samples: Sequence[Mapping[str, np.ndarray]], h: float | Sequence[float],
                  *, minimum: int = 3,
                  region_h: Mapping[str, float] | None = None,
                  zs: Sequence[float] | None = None,
                  chords: Sequence[float] | None = None,
                  log: list | None = None) -> list[ShellTopologyTemplate]:
    """``templates`` (in span order) with each region's node count set from its length.

    ``samples`` are the station polylines sampled with ``templates`` (any
    count: only their lengths are used). ``h`` is the target element edge, or
    one edge per station when the target follows the local chord. ``region_h``
    overrides it for named regions (``{"le_su": 0.02, "le_pr": 0.02}`` refines
    the leading edge). A skin region with ``chords`` given uses the smaller of
    that ceiling and the sagitta of its own polyline. Webs stay on the ceiling.
    Every region has at least ``minimum`` nodes. Stations where a region exists
    without a gap share one edge count, the most that run asks for. A shorter
    station keeps it and is logged ``block``; a station that asked for it is
    logged ``length``. A run that opens or closes takes fewer edges when the
    collapsed face cannot carry the largest ask. A section with no webs widens
    its longest region until the loop can close a rib, logged ``widened``.
    ``log`` receives one :class:`CountDecision` per region per station.
    """
    sizes = {str(k): float(v) for k, v in (region_h or {}).items()}
    if any(v <= 0.0 for v in sizes.values()):
        raise ValueError("element size h must be positive")
    if len(templates) != len(samples):
        raise ValueError("one sample set per template")
    if chords is not None and len(chords) != len(templates):
        raise ValueError("one chord per station")
    if zs is not None and len(zs) != len(templates):
        raise ValueError("one z per station")
    h_at = _station_sizes(h, len(templates))
    asked: list[dict[str, int]] = []
    widened: list[set[str]] = []
    targets_at: list[dict[str, tuple[float, str]]] = []
    lengths_at: list[dict[str, float]] = []
    for index, (tmpl, smp, hz) in enumerate(zip(templates, samples, h_at)):
        webs = {w.name for w in tmpl.webs}
        chord = None if chords is None else float(chords[index])
        taken: dict[str, int] = {}
        targets: dict[str, tuple[float, str]] = {}
        lengths: dict[str, float] = {}
        for name, pts in smp.items():
            edge, src = _region_target(name, pts, hz, chord, sizes, webs)
            targets[name] = (edge, src)
            lengths[name] = _length(pts)
            taken[name] = _even_edges(lengths[name], edge, minimum)
        grown: set[str] = set()
        if not tmpl.webs:
            before = dict(taken)
            _widen_single_cell(tmpl, smp, taken, minimum)
            grown = {name for name, count in taken.items() if count != before.get(name)}
        asked.append(taken)
        widened.append(grown)
        targets_at.append(targets)
        lengths_at.append(lengths)
    z_at = [float(t.z if zs is None else zs[i]) for i, t in enumerate(templates)]
    skin = {b.name for tmpl in templates for b in tmpl.loop}
    blocked = _block_counts(asked, lengths_at, z_at, skin, h_at)
    out: list[ShellTopologyTemplate] = []
    previous: dict[str, int] = {}
    for index, (tmpl, smp, hz) in enumerate(zip(templates, samples, h_at)):
        present = set(smp)
        for name in [name for name in previous if name not in present]:
            del previous[name]
        taken = blocked[index]
        z = float(tmpl.z if zs is None else zs[index])
        if log is not None:
            for name, pts in smp.items():
                edge, src = targets_at[index][name]
                count = int(taken[name])
                local = int(asked[index][name])
                if count != local:
                    reason = "block"
                elif name in widened[index]:
                    reason = "widened"
                else:
                    reason = "length"
                log.append(CountDecision(
                    z=z, region=name, length=_length(pts), edge=edge, edge_from=src,
                    previous=-1 if name not in previous else int(previous[name]),
                    count=count, reason=reason,
                ))
        for name in present:
            previous[name] = int(taken[name])
        least = _even_edges(0.0, hz, minimum)
        loop = tuple(replace(b, n=taken.get(b.name, least) + 1, sized=True) for b in tmpl.loop)
        webs_out = tuple(replace(w, n=taken.get(w.name, least) + 1, sized=True) for w in tmpl.webs)
        out.append(replace(tmpl, loop=loop, webs=webs_out))
    return out


__all__ = ["MIN_CELL_NODES", "CountDecision", "region_counts"]
