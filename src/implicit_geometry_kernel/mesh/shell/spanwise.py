"""Leading- and trailing-edge span lines and the sagitta checks that refine stations.

Part of :mod:`implicit_geometry_kernel.mesh.shell.shell_mesh`, split out of it in 0.9.4.
"""
from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

from implicit_geometry_kernel.blade.topology.sizing import EPS_OVER_C
from implicit_geometry_kernel.mesh.tables import (
    StationDecision,
)
from implicit_geometry_kernel.mesh.shell.stations import (
    _target_edge,
    _unique_z,
)


def _closest_on_loop_many(loop: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Closest point of a closed polyline to each row of ``pts``."""
    query = np.asarray(pts, float).reshape(-1, 3)
    ring = np.asarray(loop, float)
    if len(ring) < 2 or len(query) == 0:
        return query.copy()
    closed = np.vstack([ring, ring[:1]])
    a = closed[:-1]
    ab = np.diff(closed, axis=0)
    denom = np.maximum((ab * ab).sum(axis=1), 1e-30)
    offset = query[:, None, :] - a[None, :, :]
    t = np.clip(np.sum(offset * ab[None, :, :], axis=2) / denom, 0.0, 1.0)
    proj = a[None, :, :] + t[:, :, None] * ab[None, :, :]
    pick = np.argmin(np.linalg.norm(proj - query[:, None, :], axis=2), axis=1)
    return proj[np.arange(len(query)), pick]


def _closest_on_loop(loop: np.ndarray, xyz: np.ndarray) -> np.ndarray:
    """Closest point of a closed polyline to ``xyz``."""
    return _closest_on_loop_many(loop, np.asarray(xyz, float).reshape(1, 3))[0]


def _le_te_points(geo, zs_in: Sequence[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Placed leading and trailing edges at each ``z``, and the local chord."""
    zs = np.asarray([float(z) for z in zs_in], float)
    empty = np.zeros((0, 3), float)
    if len(zs) == 0:
        return empty, empty, np.zeros(0, float)
    chord, twist, pitch, offset, x, upper, lower = geo.interpolate_many(zs)
    signed = np.array([geo._placement_twist(float(t)) for t in twist], float)
    i_le = int(np.argmin(x))
    i_te = int(np.argmax(x))

    def at(index: int) -> np.ndarray:
        y = 0.5 * (upper[:, index] + lower[:, index])
        local_x = float(x[index]) * chord - pitch * chord
        local_y = y * chord
        angle = np.radians(signed)
        c, s = np.cos(angle), np.sin(angle)
        return np.column_stack([
            c * local_x + s * local_y + offset[:, 0],
            -s * local_x + c * local_y + offset[:, 1],
            zs,
        ])

    return at(i_le), at(i_te), np.asarray(chord, float)


def _fill_le_te_cache(geo, zs: Iterable[float], cache: dict) -> None:
    """Store leading edge, trailing edge, and chord for each missing ``z``."""
    need = [float(z) for z in zs if round(float(z), 9) not in cache]
    if not need:
        return
    le, te, chord = _le_te_points(geo, need)
    for i, z in enumerate(need):
        cache[round(float(z), 9)] = (le[i], te[i], float(chord[i]))


def _span_edges(geo, zs_in: Sequence[float], *, cache: dict | None = None):
    """Placed leading and trailing edges, and each bay's sagitta and budget.

    The edge point is the unit-chord extreme of the interpolated airfoil, placed
    with the station's chord, twist and offset. That is the point ``_le_te_point``
    used to build one station at a time. Sagitta is the bow of the straight
    ruling, so it does not build a section. ``cache`` keeps evaluated edges
    across freeze rounds.
    """
    zs = np.asarray([float(z) for z in zs_in], float)
    empty = np.zeros((0, 3), float)
    if len(zs) == 0:
        return empty, empty, np.zeros(0), np.zeros(0), np.zeros(0)
    mids = 0.5 * (zs[:-1] + zs[1:]) if len(zs) > 1 else np.zeros(0, float)
    store = cache if cache is not None else {}
    query = np.concatenate([zs, mids]) if len(mids) else zs
    _fill_le_te_cache(geo, query, store)
    le = np.stack([store[round(float(z), 9)][0] for z in zs])
    te = np.stack([store[round(float(z), 9)][1] for z in zs])
    n = len(zs)
    if n < 2:
        return le, te, np.zeros(0), np.zeros(0), np.zeros(0)
    mid_le = np.stack([store[round(float(z), 9)][0] for z in mids])
    mid_te = np.stack([store[round(float(z), 9)][1] for z in mids])
    budget = EPS_OVER_C * np.asarray([store[round(float(z), 9)][2] for z in mids], float)

    def sag(a: np.ndarray, mid: np.ndarray, b: np.ndarray) -> np.ndarray:
        ab = b - a
        length = np.linalg.norm(ab, axis=1)
        out = np.zeros(n - 1, float)
        ok = length >= 1e-12
        if np.any(ok):
            step = np.sum((mid[ok] - a[ok]) * ab[ok], axis=1) / (length[ok] * length[ok])
            out[ok] = np.linalg.norm(mid[ok] - (a[ok] + step[:, None] * ab[ok]), axis=1)
        return out

    return le, te, sag(le[:-1], mid_le, le[1:]), sag(te[:-1], mid_te, te[1:]), budget


def _le_te_point(geo, z: float, which: str) -> np.ndarray:
    """Leading or trailing edge of the shell placement at ``z``, from the loft interpolation."""
    le, te, *_rest = _span_edges(geo, (float(z),))
    return le[0] if which == "le" else te[0]


def _span_sagitta(geo, za: float, zb: float) -> tuple[float, float, float]:
    """``(le, te, budget)`` sagitta of the straight ruling against the section at mid-span."""
    _le, _te, le, te, budget = _span_edges(geo, (float(za), float(zb)))
    return float(le[0]), float(te[0]), float(budget[0])


def _sagitta_pending(geo, zs: Sequence[float], h: float | None, fraction: float | None,
                     present: set[float], *, cache: dict | None = None
                     ) -> tuple[list[tuple[float, str]], list[StationDecision]]:
    """Mid-bay stations whose leading or trailing edge bows past the budget.

    One evaluation covers every bay. A bay shorter than half the local target,
    or a mid-point that is already a station, is refused and not inserted.
    """
    _le, _te, le, te, budget = _span_edges(geo, zs, cache=cache)
    pending: list[tuple[float, str]] = []
    refused: list[StationDecision] = []
    for i in range(len(zs) - 1):
        if max(float(le[i]), float(te[i])) <= float(budget[i]):
            continue
        mid = 0.5 * (float(zs[i]) + float(zs[i + 1]))
        bay = float(zs[i + 1]) - float(zs[i])
        target = _target_edge(geo, float(zs[i + 1]), h, fraction)
        if bay < 0.5 * target:
            refused.append(StationDecision(mid, "sagitta", True, "short bay before an apex"))
            continue
        key = round(mid, 9)
        if key in present or any(abs(mid - float(s)) <= 1e-6 for s in zs):
            refused.append(StationDecision(mid, "sagitta", True, "step would land on a station"))
        else:
            pending.append((mid, "sagitta"))
            present.add(key)
    return pending, refused


def _freeze_sagitta(geo, zs: Sequence[float], h: float | None,
                    fraction: float | None) -> tuple[tuple[float, ...], dict[float, str]]:
    """Insert mid-bay stations until the leading and trailing edges hold the budget.

    This runs before any section is built. Eight rounds is the same depth the
    mesh used when it rebuilt the whole station list each pass. Edge points
    already evaluated stay in the cache; only new midpoints are interpolated.
    """
    chosen = [float(z) for z in zs]
    inserted: dict[float, str] = {}
    cache: dict = {}
    for _round in range(8):
        present = {round(z, 9) for z in chosen}
        pending, _refused = _sagitta_pending(geo, chosen, h, fraction, present, cache=cache)
        if not pending:
            break
        for z, origin in pending:
            inserted[round(float(z), 9)] = origin
        chosen = list(_unique_z(chosen + [z for z, _origin in pending]))
    return tuple(chosen), inserted


def _web_fraction_pairs(geo, zs: Sequence[float]) -> np.ndarray:
    """``(n, w, 2)`` suction and pressure chord fractions. Empty when the layout has no webs."""
    rows: list[tuple[tuple[float, float], ...]] = []
    width = 0
    for z in zs:
        feet = geo._web_feet_at(float(z))
        if feet:
            pairs = tuple((float(a), float(b)) for a, b in feet)
        else:
            layout = geo._layout_at(float(z))
            fracs = tuple(float(f) for f in getattr(layout, "web_chord_fracs", ()) or ())
            pairs = tuple((f, f) for f in fracs)
        width = max(width, len(pairs))
        rows.append(pairs)
    out = np.zeros((len(rows), width, 2), float)
    for i, pairs in enumerate(rows):
        if not pairs:
            continue
        out[i, :len(pairs)] = pairs
        if len(pairs) < width:
            out[i, len(pairs):] = pairs[-1]
    return out
