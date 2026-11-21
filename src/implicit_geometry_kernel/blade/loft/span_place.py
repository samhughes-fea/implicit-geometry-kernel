"""Placed outer-mould rings, a vertex-normal skin offset, and chord splits.

One array covers every station. A mitre join, a collapsed offset, or a cap that
changes the wall thickness stays with the section engine.
"""
from __future__ import annotations

import numpy as np


def place_oml(chord: np.ndarray, twist_deg: np.ndarray, pitch: np.ndarray, offset: np.ndarray,
              x: np.ndarray, upper: np.ndarray, lower: np.ndarray, zs: np.ndarray) -> np.ndarray:
    """Closed outer mould, one row per station. Same construction as ``oml_polyline``."""
    chord = np.asarray(chord, float).reshape(-1)
    n = len(chord)
    if n == 0:
        return np.zeros((0, 0, 3), float)
    naca_x = np.concatenate([x[::-1], x[1:]])
    naca_y = np.concatenate([np.asarray(upper, float)[:, ::-1], np.asarray(lower, float)[:, 1:]], axis=1)
    naca_x = naca_x * chord[:, None]
    naca_y = naca_y * chord[:, None]
    # ``oml_polyline`` drops a repeated trailing-edge sample before the placement.
    same = np.array([
        np.allclose(np.array([naca_x[i, 0], naca_y[i, 0]]), np.array([naca_x[i, -1], naca_y[i, -1]]))
        for i in range(n)
    ])
    if np.all(same):
        naca_x = naca_x[:, :-1]
        naca_y = naca_y[:, :-1]
    local_x = naca_x - np.asarray(pitch, float).reshape(-1)[:, None] * chord[:, None]
    local_y = naca_y
    t = np.radians(np.asarray(twist_deg, float).reshape(-1))
    c, s = np.cos(t), np.sin(t)
    off = np.asarray(offset, float).reshape(n, 2)
    xy = np.stack([
        c[:, None] * local_x + s[:, None] * local_y + off[:, 0:1],
        -s[:, None] * local_x + c[:, None] * local_y + off[:, 1:2],
    ], axis=-1)
    zcol = np.asarray(zs, float).reshape(-1)[:, None, None]
    open_ring = np.concatenate([xy, np.broadcast_to(zcol, (n, xy.shape[1], 1))], axis=-1)
    if not np.all(same) and np.any(same):
        kept = []
        for i in range(n):
            row = open_ring[i, :-1] if same[i] else open_ring[i]
            kept.append(row)
        width = max(len(row) for row in kept)
        stacked = np.empty((n, width, 3), float)
        for i, row in enumerate(kept):
            stacked[i, :len(row)] = row
            if len(row) < width:
                stacked[i, len(row):] = row[-1]
        open_ring = stacked
    return np.concatenate([open_ring, open_ring[:, :1, :]], axis=1)


def offset_rings(rings: np.ndarray, distance) -> np.ndarray:
    """Move each closed ring along its vertex normals by ``distance``.

    ``rings`` is ``(n, p, 3)``. A positive distance points toward the ring
    centroid (inward), which is half a skin thickness for the midline. Span
    ``z`` is unchanged. ``distance`` is one value or one value per station.
    """
    rings = np.asarray(rings, float)
    if rings.size == 0:
        return np.array(rings, dtype=float, copy=True)
    n = rings.shape[0]
    dist = np.asarray(distance, float)
    if dist.ndim == 0:
        dist = np.full(n, float(dist))
    else:
        dist = dist.reshape(n)
    if np.all(dist == 0.0):
        return np.array(rings, dtype=float, copy=True)
    closed = bool(np.allclose(rings[:, 0, :], rings[:, -1, :]))
    work = rings[:, :-1, :] if closed else rings
    xy = work[..., :2]
    prev = np.roll(xy, 1, axis=1)
    nxt = np.roll(xy, -1, axis=1)

    def _left(edge: np.ndarray) -> np.ndarray:
        normal = np.stack([-edge[..., 1], edge[..., 0]], axis=-1)
        length = np.linalg.norm(normal, axis=-1, keepdims=True)
        return normal / np.maximum(length, 1e-30)

    normal = _left(xy - prev) + _left(nxt - xy)
    length = np.linalg.norm(normal, axis=-1, keepdims=True)
    normal = normal / np.maximum(length, 1e-30)
    inward = xy.mean(axis=1, keepdims=True) - xy
    sign = np.sign(np.sum(normal * inward, axis=-1, keepdims=True))
    sign = np.where(sign == 0.0, 1.0, sign)
    moved = np.array(work, dtype=float, copy=True)
    moved[..., :2] = xy + dist[:, None, None] * normal * sign
    if closed:
        moved = np.concatenate([moved, moved[:, :1, :]], axis=1)
    return moved


def _interp_surface(x: np.ndarray, surface: np.ndarray, fractions: np.ndarray) -> np.ndarray:
    """``surface`` is ``(n, m)`` on a shared ``x``. ``fractions`` is ``(n, w)``."""
    x = np.asarray(x, float)
    fractions = np.clip(fractions, x[0], x[-1])
    idx = np.searchsorted(x, fractions, side="right") - 1
    idx = np.clip(idx, 0, len(x) - 2)
    x0 = x[idx]
    x1 = x[idx + 1]
    u = (fractions - x0) / np.maximum(x1 - x0, 1e-30)
    station = np.arange(fractions.shape[0])[:, None]
    return (1.0 - u) * surface[station, idx] + u * surface[station, idx + 1]


def web_lines(x: np.ndarray, upper: np.ndarray, lower: np.ndarray, chord: np.ndarray,
              pitch: np.ndarray, twist_deg: np.ndarray, offset: np.ndarray, zs: np.ndarray,
              fractions: np.ndarray) -> np.ndarray:
    """Web chords from suction to pressure at each station.

    ``fractions`` is ``(w,)`` or ``(n, w)`` when both feet share a chord
    fraction, or ``(n, w, 2)`` as ``(suction, pressure)``. The result is
    ``(n, w, 2, 3)``: the suction point and the pressure point, placed like
    the outer mould.
    """
    chord = np.asarray(chord, float).reshape(-1)
    n = len(chord)
    fractions = np.asarray(fractions, float)
    if fractions.ndim == 1:
        fractions = np.broadcast_to(fractions, (n, len(fractions), 2)).copy()
    elif fractions.ndim == 2:
        fractions = np.broadcast_to(fractions[..., None], (*fractions.shape, 2)).copy()
    elif n == 0:
        w = fractions.shape[-2] if fractions.ndim == 3 else 0
        return np.zeros((0, w, 2, 3), float)
    w = fractions.shape[1]
    if n == 0 or w == 0:
        return np.zeros((n, w, 2, 3), float)
    upper = np.asarray(upper, float)
    lower = np.asarray(lower, float)
    y_u = _interp_surface(x, upper, fractions[:, :, 0])
    y_l = _interp_surface(x, lower, fractions[:, :, 1])
    pitch = np.asarray(pitch, float).reshape(-1)
    local_x_u = fractions[:, :, 0] * chord[:, None] - pitch[:, None] * chord[:, None]
    local_x_l = fractions[:, :, 1] * chord[:, None] - pitch[:, None] * chord[:, None]
    t = np.radians(np.asarray(twist_deg, float).reshape(-1))
    c, s = np.cos(t), np.sin(t)
    off = np.asarray(offset, float).reshape(n, 2)
    z = np.asarray(zs, float).reshape(-1)

    def place(local_x: np.ndarray, local_y: np.ndarray) -> np.ndarray:
        y = local_y * chord[:, None]
        return np.stack([
            c[:, None] * local_x + s[:, None] * y + off[:, 0:1],
            -s[:, None] * local_x + c[:, None] * y + off[:, 1:2],
            np.broadcast_to(z[:, None], local_x.shape),
        ], axis=-1)

    return np.stack([place(local_x_u, y_u), place(local_x_l, y_l)], axis=2)


def split_at_arc(ring: np.ndarray, fractions: np.ndarray) -> np.ndarray:
    """Points at arc-length fractions of each closed ring.

    ``ring`` is ``(n, p, 3)`` and ``fractions`` is ``(k,)`` or ``(n, k)``.
    The result is ``(n, k, 3)``. Stations that share a node count pass the
    same fractions, so the resample is one interpolation for that group.
    """
    ring = np.asarray(ring, float)
    fractions = np.asarray(fractions, float)
    if ring.ndim != 3:
        raise ValueError("ring must be (n, p, 3)")
    n = ring.shape[0]
    if fractions.ndim == 1:
        fractions = np.broadcast_to(fractions, (n, len(fractions))).copy()
    fractions = np.clip(fractions, 0.0, 1.0)
    k = fractions.shape[1]
    if n == 0 or k == 0:
        return np.zeros((n, k, 3), float)
    closed = ring if np.allclose(ring[:, 0], ring[:, -1]) else np.concatenate([ring, ring[:, :1]], axis=1)
    step = np.linalg.norm(np.diff(closed, axis=1), axis=-1)
    length = np.concatenate([np.zeros((n, 1)), np.cumsum(step, axis=1)], axis=1)
    total = np.maximum(length[:, -1:], 1e-30)
    target = fractions * total
    # One search over every station: shift each ring onto its own span of the number line.
    span = float(np.max(total)) + 1.0
    shift = np.arange(n)[:, None] * span
    flat = np.searchsorted((length + shift).ravel(), (target + shift).ravel(), side="left") - 1
    idx = np.clip(flat.reshape(n, k) - np.arange(n)[:, None] * length.shape[1], 0, length.shape[1] - 2)
    t0 = np.take_along_axis(length, idx, axis=1)
    t1 = np.take_along_axis(length, idx + 1, axis=1)
    u = (target - t0) / np.maximum(t1 - t0, 1e-30)
    out = np.empty((n, k, 3), float)
    for axis in range(3):
        column = closed[:, :, axis]
        p0 = np.take_along_axis(column, idx, axis=1)
        p1 = np.take_along_axis(column, idx + 1, axis=1)
        out[:, :, axis] = p0 + u * (p1 - p0)
    return out


def to_chord_normal(xy: np.ndarray, chord: float, twist_deg: float, pitch: float, offset) -> np.ndarray:
    """Undo the outer placement: un-offset, un-twist, divide by chord.

    ``xy`` is ``(..., 2)`` in the blade frame. The result is the same shape,
    in chord fractions, which is the frame ``place()`` has after the pitch-axis
    map and before twist and offset.
    """
    xy = np.asarray(xy, float)
    chord = float(chord) if float(chord) else 1.0
    dx = xy[..., 0] - float(offset[0])
    dy = xy[..., 1] - float(offset[1])
    angle = np.radians(float(twist_deg))
    c, s = np.cos(angle), np.sin(angle)
    local_x = c * dx - s * dy
    local_y = s * dx + c * dy
    return np.stack([(local_x + float(pitch) * chord) / chord, local_y / chord], axis=-1)


def place_chord_normal(unit: np.ndarray, chord, twist_deg, pitch, offset, zs) -> np.ndarray:
    """Place chord-normal coordinates the way ``section_at`` places a midline.

    ``unit`` is ``(p, 2)`` or ``(n, p, 2)``. Scalars broadcast across one
    station; arrays of length ``n`` place every station. The result is
    ``(n, p, 3)``.
    """
    unit = np.asarray(unit, float)

    def _batch(value) -> np.ndarray:
        return np.atleast_1d(np.asarray(value, float)).reshape(-1)

    n = unit.shape[0] if unit.ndim == 3 else max(_batch(chord).size, _batch(zs).size, 1)
    if unit.ndim == 2:
        unit = np.broadcast_to(unit, (n, *unit.shape)).copy()
    chord = np.broadcast_to(_batch(chord), (n,)).copy()
    twist_deg = np.broadcast_to(_batch(twist_deg), (n,)).copy()
    pitch = np.broadcast_to(_batch(pitch), (n,)).copy()
    offset = np.asarray(offset, float)
    if offset.ndim == 1:
        offset = np.broadcast_to(offset, (n, 2)).copy()
    zs = np.broadcast_to(_batch(zs), (n,)).copy()
    local_x = unit[..., 0] * chord[:, None] - pitch[:, None] * chord[:, None]
    local_y = unit[..., 1] * chord[:, None]
    angle = np.radians(twist_deg)
    c, s = np.cos(angle), np.sin(angle)
    xy = np.stack([
        c[:, None] * local_x + s[:, None] * local_y + offset[:, 0:1],
        -s[:, None] * local_x + c[:, None] * local_y + offset[:, 1:2],
        np.broadcast_to(zs[:, None], local_x.shape),
    ], axis=-1)
    return xy


def resample_open_many(polylines: np.ndarray, count: int, bunch: bool = False) -> np.ndarray:
    """Arc-length resample of open polylines that share a vertex count.

    ``polylines`` is ``(n, p, d)`` and ``count`` is the sample count. One
    search covers the stack. ``bunch`` clusters samples toward the last
    vertex. Endpoints are kept. A zero-length row repeats its first point.
    """
    polylines = np.asarray(polylines, float)
    count = max(int(count), 2)
    if polylines.ndim != 3:
        raise ValueError("polylines must be (n, p, d)")
    n, p, d = polylines.shape
    out = np.empty((n, count, d), float)
    if n == 0:
        return out
    if p < 2:
        return np.repeat(polylines[:, :1, :], count, axis=1)
    step = np.linalg.norm(np.diff(polylines, axis=1), axis=-1)
    length = np.concatenate([np.zeros((n, 1)), np.cumsum(step, axis=1)], axis=1)
    total = length[:, -1]
    if bunch:
        fraction = np.sin(np.linspace(0.0, 0.5 * np.pi, count)) ** 2
    else:
        fraction = np.linspace(0.0, 1.0, count)
    target = fraction * total[:, None]
    span = float(np.max(np.abs(total))) + 1.0
    shift = np.arange(n)[:, None] * span
    flat = np.searchsorted((length + shift).ravel(), (target + shift).ravel(), side="right") - 1
    idx = flat.reshape(n, count) - np.arange(n)[:, None] * length.shape[1]
    idx = np.clip(idx, 0, p - 2)
    t0 = np.take_along_axis(length, idx, axis=1)
    t1 = np.take_along_axis(length, idx + 1, axis=1)
    span_seg = t1 - t0
    if bunch:
        denom = np.where(span_seg > 1e-15, span_seg, 1.0)
    else:
        denom = np.maximum(span_seg, 1e-15)
    u = (target - t0) / denom
    for axis in range(d):
        column = polylines[:, :, axis]
        p0 = np.take_along_axis(column, idx, axis=1)
        p1 = np.take_along_axis(column, idx + 1, axis=1)
        out[:, :, axis] = p0 + u * (p1 - p0)
    out[:, 0, :] = polylines[:, 0, :]
    out[:, -1, :] = polylines[:, -1, :]
    dead = total <= 1e-15
    if np.any(dead):
        out[dead] = polylines[dead, :1, :]
    return out
