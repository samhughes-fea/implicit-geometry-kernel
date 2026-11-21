"""Span stations of the shell mesh: target edge length, anchor and drop stations, grading.

Part of :mod:`implicit_geometry_kernel.mesh.shell.shell_mesh`, split out of it in 0.9.4.
"""
from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

from implicit_geometry_kernel.blade.topology import (
    section_cut_near,
)
from implicit_geometry_kernel.blade.topology.crossings import key_crossing_stations
from implicit_geometry_kernel.mesh.tables import (
    StationDecision,
)


TIP_H_OVER_C = 0.05
"""Default shell edge as a fraction of the local chord. ``h`` or ``n_span`` leaves it off."""


def _unique_z(zs: Iterable[float], *, abs_tol: float = 1e-3) -> tuple[float, ...]:
    """Sorted stations merged when closer than ``abs_tol`` (metres).

    A key-crossing can land within a millimetre of an anchor; those must
    collapse to one station so neighbouring bays weld the same nodes.
    """
    ordered = sorted(float(z) for z in zs)
    if not ordered:
        return ()
    out = [ordered[0]]
    for z in ordered[1:]:
        if abs(z - out[-1]) > abs_tol:
            out.append(z)
    return tuple(out)


def _chord_at(geo, z: float) -> float:
    return float(geo._span_scalar(getattr(geo, "_chord", None),
                                  [s.chord for s in geo.stations], float(z)))


def _target_edge(geo, z: float, h: float | None, h_over_c: float | None) -> float:
    """Local quad edge: ``h`` in metres, ``h_over_c`` times the local chord, or the smaller of the two."""
    s = None
    if h_over_c is not None:
        if float(h_over_c) <= 0.0:
            raise ValueError("h_over_c must be positive")
        s = float(h_over_c) * _chord_at(geo, z)
    if h is not None:
        if float(h) <= 0.0:
            raise ValueError("element size h must be positive")
        s = float(h) if s is None else min(float(h), s)
    if s is None or s <= 0.0:
        raise ValueError("element size must be positive")
    return s


def _n_bays(distance: float, s: float) -> int:
    """Equal bays whose step is closest to ``s`` and still inside ``[s/2, s]``.

    A gap of at most ``s`` is one bay, including a short gap between two anchors.
    """
    if s <= 0.0 or distance <= s * (1.0 + 1e-9):
        return 1
    lo = max(int(np.ceil(distance / s - 1e-9)), 1)
    hi = int(np.floor(2.0 * distance / s + 1e-9))
    if hi < lo:
        return lo
    target = distance / s
    return min(range(lo, hi + 1), key=lambda n: (abs(n - target), n))


def _fill_gaps(stations: Sequence[float], s_at) -> tuple[float, ...]:
    """Fit regular stations into the gaps. ``s_at(z)`` is the local target edge."""
    zs = sorted(float(z) for z in stations)
    if len(zs) < 2:
        return tuple(zs)
    out = [zs[0]]
    for z in zs[1:]:
        start = out[-1]
        s = min(s_at(start), s_at(z))
        n = _n_bays(z - start, s)
        if n > 1:
            step = (z - start) / n
            out.extend(start + step * k for k in range(1, n))
        out.append(z)
    return _unique_z(out)


def _drop_stations(geo) -> tuple[float, ...]:
    """Span where a web with no stored end leaves the section.

    Leading-edge and pitch-axis webs are dropped, not slid, when they miss the
    inner skin. The station sits on that change.
    """
    webs = [w for w in getattr(geo.region_graph, "webs", ()) or () if w.z1 is None]
    if not webs:
        return ()
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    tip = section_cut_near(geo, z1)
    if tip is None:
        return ()
    dropped_tip = {int(i) for i in (getattr(tip, "dropped_web_indices", ()) or ())}
    need = [w for w in webs if int(w.index) in dropped_tip]
    if not need:
        return ()
    root = section_cut_near(geo, z0)
    dropped_root = {int(i) for i in (getattr(root, "dropped_web_indices", ()) or ())} if root is not None else set()
    found: list[float] = []
    span = max(z1 - z0, 1.0)
    for web in need:
        if int(web.index) in dropped_root:
            continue
        lo, hi = z0, z1
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            cut = section_cut_near(geo, mid)
            dropped = {int(i) for i in (getattr(cut, "dropped_web_indices", ()) or ())} if cut is not None else set()
            if int(web.index) in dropped:
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-6 * span:
                break
        if z0 + 1e-9 < hi < z1 - 1e-9:
            found.append(hi)
    return tuple(found)


def _anchor_stations(geo) -> tuple[float, ...]:
    """Spanwise stations the section sample must contain.

    The root, the tip and every defining station. The birth and death of each
    web, including a web that leaves the section with no authored end. Every
    rib and other transverse member, the ends of stringers and ply drops, and
    a key crossing. The section sample is built from this list, and the mesh
    uses the same stations as its anchors.
    """
    cached = geo.__dict__.get("_mesh_anchor_z")
    if cached is not None:
        return cached
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = [float(station.z) for station in geo.stations]
    zs.append(float(geo.tip_closure_z))
    zs.extend(float(z) for z in geo.forced_stations())
    zs.extend(float(member.z) for member in geo.transverse)
    for anchor in getattr(geo.region_graph, "anchors", ()) or ():
        for z in (getattr(anchor, "z0", None), getattr(anchor, "z1", None)):
            if z is not None and z0 < float(z) < z1:
                zs.append(float(z))
    zs.extend(_drop_stations(geo))
    found = _unique_z(z for z in zs if z0 - 1e-9 <= z <= z1 + 1e-9)
    found = _unique_z(list(found) + list(key_crossing_stations(geo, found)))
    geo.__dict__["_mesh_anchor_z"] = found
    return found


def _refine_stations(zs: Sequence[float], z0: float, z1: float,
                     refine: Sequence[tuple[float, float, float]]) -> tuple[float, ...]:
    out = set(zs)
    for a, b, dz in refine:
        a, b = max(float(min(a, b)), z0), min(float(max(a, b)), z1)
        if float(dz) <= 0.0:
            raise ValueError("refinement spacing dz must be positive")
        if b > a:
            count = int(np.ceil((b - a) / float(dz) - 1e-9)) + 1
            out.update(float(z) for z in np.linspace(a, b, count))
    return _unique_z(out)


def _median_skin_edge(tmpl, smp) -> float | None:
    """Median chordwise edge of the skin regions. Webs and the blunt wrap are left out."""
    lengths = []
    for branch in tmpl.loop:
        if branch.wrap or branch.name not in smp:
            continue
        pts = np.asarray(smp[branch.name], float)
        edges = len(pts) - 1
        if edges <= 0:
            continue
        length = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
        if length > 0.0:
            lengths.append(length / edges)
    if not lengths:
        return None
    return float(np.median(lengths))


def _grade_stations(geo, zs: Sequence[float], tmpls, smps, h: float | None,
                    h_over_c: float | None):
    """Split a bay longer than about 1.25 times the outboard skin edge. One pass.

    The new step stays at least half the local target, so a station is not
    placed on top of an anchor. A skin edge already shorter than that half is
    the geometry (a narrow cap, a short web), and the span step does not chase it.
    A bay left long is returned as a refused grade, with the reason.
    """
    anchors = {round(z, 9) for z in _anchor_stations(geo)}
    out = [float(zs[0])]
    changed = False
    refused: list[StationDecision] = []
    for i in range(len(zs) - 1):
        za, zb = float(zs[i]), float(zs[i + 1])
        med = _median_skin_edge(tmpls[i + 1], smps[i + 1])
        bay = zb - za
        mid = 0.5 * (za + zb)
        if med is not None and bay > 1.25 * med * (1.0 + 1e-12):
            step_min = 0.5 * _target_edge(geo, zb, h, h_over_c)
            if med < step_min:
                refused.append(StationDecision(
                    mid, "grade", True, "skin edge already under half the target"))
            else:
                n = max(int(np.ceil(bay / med - 1e-9)), 1)
                step = bay / n
                if n > 1 and step + 1e-9 >= step_min:
                    clear = True
                    for k in range(1, n):
                        p = za + step * k
                        for end in (za, zb):
                            if round(end, 9) in anchors and abs(p - end) < step_min - 1e-8:
                                clear = False
                    if clear:
                        out.extend(za + step * k for k in range(1, n))
                        changed = True
                    else:
                        refused.append(StationDecision(
                            mid, "grade", True, "step would land on an anchor"))
                else:
                    refused.append(StationDecision(
                        mid, "grade", True, "split step under half the target"))
        out.append(zb)
    return _unique_z(out), changed, tuple(refused)


def mesh_stations(geo, n_span: int | None = None, *, h: float | None = None,
                  h_over_c: float | None = None,
                  refine: Sequence[tuple[float, float, float]] = ()) -> tuple[float, ...]:
    """Span stations of the mesh.

    With ``h`` (metres) or ``h_over_c`` (a fraction of the local chord), the
    anchors are the defining stations, the root, the tip and the topology
    events, and regular stations fill the gaps so each bay is in ``[s/2, s]``
    of the local target edge.
    With ``n_span`` and no size, the stations are that count, plus the same
    anchors. ``refine`` holds ``(z0, z1, dz)`` windows. A key crossing found
    between stations is kept, and that bay is filled again.
    """
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    if h is not None and float(h) <= 0.0:
        raise ValueError("element size h must be positive")
    if h_over_c is not None and float(h_over_c) <= 0.0:
        raise ValueError("h_over_c must be positive")
    anchors = _anchor_stations(geo)
    if h is None and h_over_c is None:
        if n_span is None:
            raise ValueError("give n_span, h, or h_over_c")
        zs = _unique_z(list(np.linspace(z0, z1, int(n_span))) + list(anchors))
        zs = _refine_stations(zs, z0, z1, refine)
        base = [z for z in zs if z0 - 1e-9 <= z <= z1 + 1e-9]
        return _unique_z(list(base) + list(key_crossing_stations(geo, base)))
    s_at = lambda z, _h=h, _f=h_over_c: _target_edge(geo, z, _h, _f)
    filled = _fill_gaps(anchors, s_at)
    filled = _refine_stations(filled, z0, z1, refine)
    crossings = key_crossing_stations(geo, filled)
    if crossings:
        filled = _fill_gaps(tuple(filled) + tuple(crossings), s_at)
    return filled


def _size_fraction(n_span: int | None, h: float | None, h_over_c: float | None,
                   explicit: bool) -> float | None:
    """The chord fraction in force, or ``None`` when the caller sized by metres or by count.

    A bare call uses :data:`TIP_H_OVER_C`. ``h`` alone, ``n_span``, and an
    explicit station list leave that default off. A passed ``h_over_c`` is the
    fraction and overrides the default.
    """
    if h_over_c is not None:
        if float(h_over_c) <= 0.0:
            raise ValueError("h_over_c must be positive")
        return float(h_over_c)
    if h is not None or n_span is not None or explicit:
        return None
    return TIP_H_OVER_C
