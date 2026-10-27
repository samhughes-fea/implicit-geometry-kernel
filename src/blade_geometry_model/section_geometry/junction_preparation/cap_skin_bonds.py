"""Bond every spar cap to the skin arcs it lies on (thin_wall_section TW-001).

A cap is bonded along its length to the skin of its own surface. The thin-wall
kernel needs that bond as a **vertex correspondence**: every cap segment lies
against exactly one skin segment, so the cap's axial rate enters the skin
linearly along each pair and the skin's shear flow stays exactly quadratic.

This pass makes the correspondence exist and emits it:

* every cap vertex is paired with its perpendicular foot on the skin half-ring
  of the cap's surface; the cap-to-skin map is piecewise linear between those
  pairs, and every skin vertex inside the cap's footprint takes its cap partner
  from that same map (one map, so the pairing is one-to-one by construction);
* a foot that lands close to an existing vertex of the other strip is paired
  with that vertex instead (``_PAIR_MERGE_FRACTION`` of the adjacent segment) —
  this chooses *which* points correspond, never where anything is; every other
  foot is inserted as a new vertex **on** its polyline, so no geometry moves;
* the correspondence must be monotone along both strips, or the cap is refused.

Which skin a cap belongs to is declared (``ShellMidlineStrip.surface`` against
``skin_side``), not inferred. The result is carried on the cap as
``ShellMidlineStrip.bonds``, one :class:`StripBond` per host skin arc.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

from blade_geometry_model.section_contract.shell_midline_strip import (
    ShellMidlineStrip,
    StripBond,
)

from ..midlines.cap_skin_ray_resample import _surface_side

__all__ = ["bond_caps_to_skin"]

_JUNCTION_TOL_M = 1e-9
"""Adjacent skin arcs share their ends to this distance (they are bit-identical)."""

_PAIR_MERGE_FRACTION = 0.25
"""A foot within this fraction of the adjacent segment of an existing vertex is
paired with that vertex rather than inserted beside it."""


@dataclass
class _Arc:
    strip_index: int
    strip: ShellMidlineStrip
    reversed: bool

    @property
    def oriented(self) -> NDArray[np.float64]:
        points = np.asarray(self.strip.midline_twisted_pitch_axis_blade_frame, dtype=np.float64)
        return points[::-1] if self.reversed else points


def _chain(arcs: list[tuple[int, ShellMidlineStrip]]) -> list[_Arc]:
    """The skin arcs of one surface, head to tail."""
    pieces = [(i, s, np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=np.float64)) for i, s in arcs]
    ends = [(p[0], p[-1]) for _i, _s, p in pieces]

    def meets(point, index):
        return any(
            np.hypot(*(point - end)) <= _JUNCTION_TOL_M
            for other, pair in enumerate(ends)
            if other != index
            for end in pair
        )

    open_ends = [(k, 0) for k in range(len(pieces)) if not meets(ends[k][0], k)]
    open_ends += [(k, 1) for k in range(len(pieces)) if not meets(ends[k][1], k)]
    if len(open_ends) != 2:
        raise ValueError(f"the skin arcs of one surface do not form a single chain ({len(open_ends)} open ends)")
    start, which = open_ends[0]
    chain = [_Arc(pieces[start][0], pieces[start][1], reversed=(which == 1))]
    used = {start}
    while len(used) < len(pieces):
        tail = chain[-1].oriented[-1]
        for k, (_i, _s, p) in enumerate(pieces):
            if k in used:
                continue
            if np.hypot(*(p[0] - tail)) <= _JUNCTION_TOL_M:
                chain.append(_Arc(pieces[k][0], pieces[k][1], reversed=False))
                break
            if np.hypot(*(p[-1] - tail)) <= _JUNCTION_TOL_M:
                chain.append(_Arc(pieces[k][0], pieces[k][1], reversed=True))
                break
        else:
            raise ValueError("the skin arcs of one surface do not chain head to tail")
        used.add(k)
    return chain


def _prefix(points: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(points, axis=0).T))])


def _foot(points: NDArray[np.float64], prefix: NDArray[np.float64], target: NDArray[np.float64]):
    """``(arc parameter, foot point, segment, t)`` of the closest point on a polyline."""
    a, b = points[:-1], points[1:]
    ab = b - a
    t = np.clip(np.einsum("ij,ij->i", target - a, ab) / np.einsum("ij,ij->i", ab, ab), 0.0, 1.0)
    feet = a + t[:, None] * ab
    segment = int(np.argmin(np.hypot(*(feet - target).T)))
    return float(prefix[segment] + t[segment] * np.hypot(*ab[segment])), feet[segment], segment, float(t[segment])


def _merge_distance(points: NDArray[np.float64], vertex: int) -> float:
    lengths = []
    if vertex > 0:
        lengths.append(np.hypot(*(points[vertex] - points[vertex - 1])))
    if vertex < len(points) - 1:
        lengths.append(np.hypot(*(points[vertex + 1] - points[vertex])))
    return _PAIR_MERGE_FRACTION * float(min(lengths))


def _with_inserts(points: NDArray[np.float64], inserts: list[tuple[int, float, NDArray[np.float64]]]):
    """``points`` with ``(segment, t, point)`` inserted in order."""
    by_segment: dict[int, list[tuple[float, NDArray[np.float64]]]] = {}
    for segment, t, point in inserts:
        by_segment.setdefault(segment, []).append((t, point))
    out = [points[0]]
    for segment in range(len(points) - 1):
        for _t, point in sorted(by_segment.get(segment, []), key=lambda item: item[0]):
            out.append(point)
        out.append(points[segment + 1])
    return np.asarray(out, dtype=np.float64)


def _bond_one_cap(cap: ShellMidlineStrip, chain: list[_Arc]):
    """New cap points, new points per chain arc, and the bonds.

    One map, so the correspondence is one-to-one by construction: each cap
    vertex is paired with its perpendicular foot on the skin chain (and the feet
    must advance monotonically), and the cap-to-skin map is piecewise linear
    between those pairs. A skin vertex inside the footprint takes its cap partner
    from that same map -- a second, independent projection back onto the cap is
    not the inverse of the first and can interleave near a junction.
    """
    cap_points = np.asarray(cap.midline_twisted_pitch_axis_blade_frame, dtype=np.float64)
    chain_points = np.vstack([chain[0].oriented] + [arc.oriented[1:] for arc in chain[1:]])
    chain_prefix = _prefix(chain_points)
    arc_bounds = np.cumsum([0] + [len(arc.oriented) - 1 for arc in chain])
    cap_prefix = _prefix(cap_points)

    feet = [_foot(chain_points, chain_prefix, point) for point in cap_points]
    sigma = np.array([f[0] for f in feet])
    steps = np.diff(sigma)
    if not (np.all(steps > 0.0) or np.all(steps < 0.0)):
        raise ValueError(f"cap {cap.label!r} does not run monotonically along its skin")
    chain_partner: list[NDArray[np.float64]] = [f[1] for f in feet]
    chain_insert: list[tuple[int, float] | None] = [(f[2], f[3]) for f in feet]
    paired_chain_vertex: dict[int, int] = {}

    low, high = float(sigma.min()), float(sigma.max())
    inside = [j for j in range(len(chain_points)) if low - 1e-12 <= chain_prefix[j] <= high + 1e-12]
    # A skin vertex beside a cap vertex's foot is paired with that cap vertex.
    for j in inside:
        i = int(np.argmin(np.abs(sigma - chain_prefix[j])))
        if i in paired_chain_vertex.values():
            continue
        if abs(sigma[i] - chain_prefix[j]) <= _merge_distance(chain_points, j):
            paired_chain_vertex[j] = i
            sigma[i] = chain_prefix[j]
            chain_partner[i] = chain_points[j]
            chain_insert[i] = None

    pairs = [
        {"sigma": float(sigma[i]), "chain": chain_partner[i], "cap": cap_points[i], "chain_insert": chain_insert[i],
         "cap_insert": None}
        for i in range(len(cap_points))
    ]
    order = np.argsort(sigma)
    sigma_sorted, tau_sorted = sigma[order], cap_prefix[order]
    for j in inside:
        if j in paired_chain_vertex or not (low < chain_prefix[j] < high):
            continue
        # The cap point the map sends this skin vertex to: interpolate arc length.
        tau = float(np.interp(chain_prefix[j], sigma_sorted, tau_sorted))
        segment = int(np.clip(np.searchsorted(cap_prefix, tau, side="right") - 1, 0, len(cap_points) - 2))
        span = cap_prefix[segment + 1] - cap_prefix[segment]
        t = (tau - cap_prefix[segment]) / span
        point = cap_points[segment] + t * (cap_points[segment + 1] - cap_points[segment])
        pairs.append({"sigma": float(chain_prefix[j]), "chain": chain_points[j], "cap": point, "chain_insert": None,
                      "cap_insert": (segment, float(t)), "tau": tau})

    pairs.sort(key=lambda p: p["sigma"])
    for index, p in enumerate(pairs):
        p.setdefault("tau", None)
    taus = [p["tau"] if p["tau"] is not None else float(np.interp(p["sigma"], sigma_sorted, tau_sorted)) for p in pairs]
    if not (np.all(np.diff(taus) > 0.0) or np.all(np.diff(taus) < 0.0)) or not np.all(
        np.diff([p["sigma"] for p in pairs]) > 0.0
    ):
        raise ValueError(f"cap {cap.label!r}: the cap-to-skin correspondence is not one-to-one")

    new_cap = _with_inserts(cap_points, [(*p["cap_insert"], p["cap"]) for p in pairs if p["cap_insert"] is not None])
    arc_inserts: dict[int, list[tuple[int, float, NDArray[np.float64]]]] = {}
    for p in pairs:
        if p["chain_insert"] is None:
            continue
        segment, t = p["chain_insert"]
        arc = int(np.searchsorted(arc_bounds, segment, side="right") - 1)
        local = segment - int(arc_bounds[arc])
        if chain[arc].reversed:
            count = len(chain[arc].oriented) - 1
            local, t = count - 1 - local, 1.0 - t
        arc_inserts.setdefault(arc, []).append((local, t, p["chain"]))

    bonds: list[StripBond] = []
    for arc_index, arc in enumerate(chain):
        lo, hi = chain_prefix[arc_bounds[arc_index]], chain_prefix[arc_bounds[arc_index + 1]]
        members = [p for p in pairs if lo - 1e-12 <= p["sigma"] <= hi + 1e-12]
        if len(members) >= 2:
            bonds.append(
                StripBond(
                    host_label=str(arc.strip.label),
                    points=np.asarray([p["cap"] for p in members], dtype=np.float64),
                    host_points=np.asarray([p["chain"] for p in members], dtype=np.float64),
                )
            )
    return new_cap, arc_inserts, tuple(bonds)


def bond_caps_to_skin(strips: tuple[ShellMidlineStrip, ...]) -> tuple[ShellMidlineStrip, ...]:
    """Pair every cap with its skin, inserting correspondence vertices; emit the bonds.

    Raises
    ------
    ValueError
        A surface's skin arcs do not chain, or a cap does not run monotonically
        along its skin.
    """
    out = list(strips)
    for cap_index, cap in enumerate(strips):
        if cap.kind != "cap":
            continue
        side = _surface_side(cap.surface)
        if side is None:
            continue
        hosts = [(i, s) for i, s in enumerate(out) if s.kind == "skin" and _surface_side(s.skin_side) == side]
        if not hosts:
            continue
        chain = _chain(hosts)
        new_cap, arc_inserts, bonds = _bond_one_cap(out[cap_index], chain)
        for arc_index, inserts in arc_inserts.items():
            arc = chain[arc_index]
            points = np.asarray(out[arc.strip_index].midline_twisted_pitch_axis_blade_frame, dtype=np.float64)
            out[arc.strip_index] = replace(
                out[arc.strip_index], midline_twisted_pitch_axis_blade_frame=_with_inserts(points, inserts)
            )
        out[cap_index] = replace(out[cap_index], midline_twisted_pitch_axis_blade_frame=new_cap, bonds=bonds)
    return tuple(out)
