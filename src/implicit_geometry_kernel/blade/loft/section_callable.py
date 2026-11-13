"""Loft an outer shape given as ``section(z)`` until a linear loft meets a tolerance.

The callable is already the consumer's blend (planform laws, a root blend). Stations
are added until the straight loft between them lies within ``tolerance`` metres of
``section(z)``, measured on the placed outer mould line. Interpolation between those
stations is linear in the sampled surfaces and in chord, twist, pitch axis and offset.
It does not run :func:`_morph_surfaces` again: that blend is for station-table blades.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence, Union

import numpy as np

from implicit_geometry_kernel.blade.airfoil import Airfoil
from implicit_geometry_kernel.blade.cut import _twist_offset

SectionCallable = Callable[[float], object]


@dataclass(frozen=True)
class SectionSample:
    """An already-evaluated outer section at one span station.

    ``x``, ``upper`` and ``lower`` are unit chord (leading edge at 0). ``twist`` is
    degrees, nose up. ``pivot_fraction`` is the pitch-axis chord fraction from the
    leading edge. ``offset`` is the ``(x, y)`` shift of that point, in metres.
    """

    z: float
    chord: float
    twist: float
    pivot_fraction: float
    offset: tuple[float, float]
    x: np.ndarray
    upper: np.ndarray
    lower: np.ndarray


def cosine_x(n: int) -> np.ndarray:
    """The chord fractions ``interpolate`` samples a surface on."""
    return 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, int(n))))


def airfoil_from_surfaces(x: np.ndarray, upper: np.ndarray, lower: np.ndarray, name: str = "section") -> Airfoil:
    """An ``Airfoil`` whose Selig contour is the given unit-chord surfaces."""
    xu = np.asarray(x, float)
    pts = np.vstack([
        np.column_stack([xu[::-1], np.asarray(upper, float)[::-1]]),
        np.column_stack([xu[1:], np.asarray(lower, float)[1:]]),
    ])
    return Airfoil(pts, name=name)


def placed_oml(sample: SectionSample) -> np.ndarray:
    """Closed outer mould line of ``sample`` in the blade frame, shape ``(N, 2)``."""
    x = np.asarray(sample.x, float)
    upper = np.asarray(sample.upper, float)
    lower = np.asarray(sample.lower, float)
    naca = np.vstack([
        np.column_stack([x[::-1], upper[::-1]]),
        np.column_stack([x[1:], lower[1:]]),
    ]) * float(sample.chord)
    if np.allclose(naca[0], naca[-1]):
        naca = naca[:-1]
    local = np.column_stack([naca[:, 0] - float(sample.pivot_fraction) * float(sample.chord), naca[:, 1]])
    xy = _twist_offset(local, float(sample.twist), sample.offset)
    if np.linalg.norm(xy[0] - xy[-1]) > 1e-12:
        xy = np.vstack([xy, xy[0]])
    return xy


def polyline_gap(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric distance between two polylines: the furthest point of either from the other."""
    aa = np.asarray(a, float).reshape(-1, 2)
    bb = np.asarray(b, float).reshape(-1, 2)
    if len(aa) < 2 or len(bb) < 2:
        raise ValueError("a polyline gap needs two polylines of at least two points")
    return float(max(_max_distance(aa, bb), _max_distance(bb, aa)))


def _max_distance(points: np.ndarray, poly: np.ndarray) -> float:
    if np.linalg.norm(poly[0] - poly[-1]) > 1e-12:
        poly = np.vstack([poly, poly[0]])
    a = poly[:-1]
    b = poly[1:]
    ab = b - a
    length2 = np.sum(ab * ab, axis=1)
    ap = points[:, None, :] - a[None, :, :]
    t = np.clip(np.sum(ap * ab, axis=2) / np.maximum(length2, 1e-30), 0.0, 1.0)
    proj = a[None, :, :] + t[:, :, None] * ab[None, :, :]
    return float(np.linalg.norm(points[:, None, :] - proj, axis=2).min(axis=1).max())


def lerp_sample(a: SectionSample, b: SectionSample, w: float, z: float) -> SectionSample:
    """Linear loft of two samples. Surfaces are not remorphed."""
    w = float(np.clip(w, 0.0, 1.0))
    offset = (
        (1.0 - w) * float(a.offset[0]) + w * float(b.offset[0]),
        (1.0 - w) * float(a.offset[1]) + w * float(b.offset[1]),
    )
    return SectionSample(
        z=float(z),
        chord=(1.0 - w) * float(a.chord) + w * float(b.chord),
        twist=(1.0 - w) * float(a.twist) + w * float(b.twist),
        pivot_fraction=(1.0 - w) * float(a.pivot_fraction) + w * float(b.pivot_fraction),
        offset=offset,
        x=np.asarray(a.x, float),
        upper=(1.0 - w) * np.asarray(a.upper, float) + w * np.asarray(b.upper, float),
        lower=(1.0 - w) * np.asarray(a.lower, float) + w * np.asarray(b.lower, float),
    )


def _as_sample(value, z: float, x: np.ndarray) -> SectionSample:
    if isinstance(value, SectionSample):
        return _resample(value, x)
    if isinstance(value, dict):
        shape = value.get("airfoil", value.get("surfaces"))
        return _from_parts(
            z,
            value["chord"],
            value.get("twist", 0.0),
            value.get("pivot_fraction", value.get("pitch_axis")),
            value.get("offset", (0.0, 0.0)),
            shape,
            x,
        )
    if isinstance(value, (tuple, list)) and len(value) == 5:
        chord, twist, pivot, offset, shape = value
        return _from_parts(z, chord, twist, pivot, offset, shape, x)
    raise TypeError(
        "section(z) must return a SectionSample, a mapping, or "
        "(chord, twist, pivot_fraction, offset, airfoil or (x, upper, lower))"
    )


def _resample(sample: SectionSample, x: np.ndarray) -> SectionSample:
    src = np.asarray(sample.x, float)
    grid = np.asarray(x, float)
    if src.shape == grid.shape and np.allclose(src, grid):
        upper = np.asarray(sample.upper, float)
        lower = np.asarray(sample.lower, float)
    else:
        upper = np.interp(grid, src, np.asarray(sample.upper, float))
        lower = np.interp(grid, src, np.asarray(sample.lower, float))
    return SectionSample(
        z=float(sample.z),
        chord=float(sample.chord),
        twist=float(sample.twist),
        pivot_fraction=float(sample.pivot_fraction),
        offset=(float(sample.offset[0]), float(sample.offset[1])),
        x=grid,
        upper=upper,
        lower=lower,
    )


def _from_parts(z, chord, twist, pivot, offset, shape, x: np.ndarray) -> SectionSample:
    if pivot is None:
        raise ValueError("section(z) needs a pivot_fraction")
    grid = np.asarray(x, float)
    if isinstance(shape, Airfoil):
        upper = np.asarray(shape.surface("upper", grid), float)
        lower = np.asarray(shape.surface("lower", grid), float)
    else:
        sx, upper, lower = shape
        return _resample(
            SectionSample(
                float(z), float(chord), float(twist), float(pivot),
                (float(offset[0]), float(offset[1])),
                np.asarray(sx, float), np.asarray(upper, float), np.asarray(lower, float),
            ),
            grid,
        )
    return SectionSample(
        z=float(z),
        chord=float(chord),
        twist=float(twist),
        pivot_fraction=float(pivot),
        offset=(float(offset[0]), float(offset[1])),
        x=grid,
        upper=upper,
        lower=lower,
    )


def refine_section_samples(
    section: SectionCallable,
    z: Union[float, Sequence[float]],
    tolerance: float,
    x: np.ndarray,
    *,
    knots: Sequence[float] = (),
    probes: Sequence[float] = (0.25, 0.5, 0.75),
    max_stations: int = 2000,
) -> tuple[SectionSample, ...]:
    """Samples of ``section`` from the first span to the last, including ``knots``.

    An interval is kept when the linear loft at each probe is within ``tolerance``
    metres of ``section`` at that span. Otherwise a station is added at the worst
    probe and the two sides are refined. The callable is evaluated once per station.
    """
    tol = float(tolerance)
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tolerance must be a positive length in metres")
    span = [float(v) for v in np.atleast_1d(np.asarray(z, float)).reshape(-1)]
    span.extend(float(v) for v in knots)
    if len(span) < 2:
        raise ValueError("a section loft needs at least two span positions")
    by_key: dict[float, float] = {}
    for value in span:
        by_key.setdefault(round(value, 9), value)
    knots_z = [by_key[key] for key in sorted(by_key)]
    if any(np.diff(knots_z) <= 0.0):
        raise ValueError("section span positions must be strictly increasing")
    cache: dict[float, SectionSample] = {}

    def at(zz: float) -> SectionSample:
        key = round(float(zz), 9)
        hit = cache.get(key)
        if hit is None:
            hit = _as_sample(section(float(zz)), float(zz), x)
            if abs(hit.z - float(zz)) > 1e-6 * max(1.0, abs(float(zz))):
                hit = SectionSample(
                    float(zz), hit.chord, hit.twist, hit.pivot_fraction, hit.offset,
                    hit.x, hit.upper, hit.lower,
                )
            cache[key] = hit
        return hit

    by_z = {round(s.z, 9): s for s in (at(zz) for zz in knots_z)}
    span_length = max(by_z) - min(by_z)
    min_dz = 1e-9 * max(1.0, abs(span_length))
    probe_w = tuple(float(p) for p in probes)
    if not probe_w or any(not 0.0 < p < 1.0 for p in probe_w):
        raise ValueError("probes must lie inside (0, 1)")
    pending = [
        (by_z[a].z, by_z[b].z)
        for a, b in zip(sorted(by_z), sorted(by_z)[1:])
    ]

    def _probe_gap(a: SectionSample, b: SectionSample, w: float) -> tuple[float, float]:
        z_probe = a.z + w * (b.z - a.z)
        gap = polyline_gap(placed_oml(lerp_sample(a, b, w, z_probe)), placed_oml(at(z_probe)))
        return gap, z_probe

    while pending:
        if len(by_z) > max_stations:
            raise ValueError(
                f"section loft exceeded {max_stations} stations before reaching {tol:g} m"
            )
        z0, z1 = pending.pop()
        a, b = by_z[round(z0, 9)], by_z[round(z1, 9)]
        if b.z <= a.z:
            raise ValueError("refined stations crossed")
        if b.z - a.z <= min_dz:
            gap, z_probe = _probe_gap(a, b, 0.5)
            if gap > tol:
                raise ValueError(
                    f"section loft at z = {z_probe:g} is {gap:.3g} m "
                    f"from the linear loft, above {tol:g} m"
                )
            continue
        worst = 0.0
        worst_z = 0.5 * (a.z + b.z)
        for w in probe_w:
            gap, z_probe = _probe_gap(a, b, w)
            if gap > worst:
                worst = gap
                worst_z = z_probe
        if worst <= tol:
            continue
        added = at(worst_z)
        if not a.z < added.z < b.z:
            added = at(0.5 * (a.z + b.z))
        if not a.z < added.z < b.z or round(added.z, 9) in (round(a.z, 9), round(b.z, 9)):
            raise ValueError(
                f"section loft at z = {worst_z:g} is {worst:.3g} m "
                f"from the linear loft, above {tol:g} m"
            )
        by_z[round(added.z, 9)] = added
        pending.append((a.z, added.z))
        pending.append((added.z, b.z))
    samples = tuple(by_z[key] for key in sorted(by_z))
    if any(np.diff([s.z for s in samples]) <= 0.0):
        raise ValueError("refined stations are not strictly increasing")
    return samples


__all__ = [
    "SectionCallable",
    "SectionSample",
    "airfoil_from_surfaces",
    "cosine_x",
    "lerp_sample",
    "placed_oml",
    "polyline_gap",
    "refine_section_samples",
]
