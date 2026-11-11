"""Stations where two skin keys swap order along the span.

A cap edge can pass a web foot between two stations (the cap narrows past a
web, or the web runs at a different chord fraction from the cap). The region
between the two keys then shrinks to nothing and grows again reversed. A
ruled strip across such a band overlaps its neighbours, so the mesh needs a
station at the crossing itself: there the region has no length, and the loft
sees one region end and another begin, which it already meshes conformingly.
"""
from __future__ import annotations

from typing import Sequence

from implicit_geometry_kernel.blade.topology.ribs import _station_side_fracs, section_cut_near

_SIDES = ("suction", "pressure")


def _side_fracs_at(geo, z: float) -> dict[str, dict[str, float]] | None:
    cut = section_cut_near(geo, float(z))
    if cut is None:
        return None
    return _station_side_fracs(geo.region_graph, cut, cut.shell_graph(), float(cut.z), True)


def _swapped_pairs(a: dict[str, dict[str, float]], b: dict[str, dict[str, float]]
                   ) -> list[tuple[str, str, str]]:
    """``(side, key_1, key_2)`` for every key pair whose order differs between ``a`` and ``b``."""
    out = []
    for side in _SIDES:
        fa, fb = a.get(side, {}), b.get(side, {})
        names = sorted(n for n in fa if n in fb)
        for i, n1 in enumerate(names):
            for n2 in names[i + 1:]:
                if (fa[n1] - fa[n2]) * (fb[n1] - fb[n2]) < 0.0:
                    out.append((side, n1, n2))
    return out


def _crossing_z(geo, side: str, n1: str, n2: str, z0: float, z1: float, d0: float,
                tol: float) -> float | None:
    """Bisect for the ``z`` where keys ``n1`` and ``n2`` meet on ``side``."""
    lo, hi, d_lo = float(z0), float(z1), float(d0)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        fr = _side_fracs_at(geo, mid)
        if fr is None or n1 not in fr.get(side, {}) or n2 not in fr.get(side, {}):
            return None
        d = fr[side][n1] - fr[side][n2]
        if d == 0.0:
            return mid
        if (d < 0.0) == (d_lo < 0.0):
            lo, d_lo = mid, d
        else:
            hi = mid
    return 0.5 * (lo + hi)


def key_crossing_stations(geo, zs: Sequence[float], *, rel_tol: float = 1e-9) -> tuple[float, ...]:
    """Spans between consecutive ``zs`` where two skin keys swap order, to ``rel_tol`` of the span.

    Every swap found between two stations adds one station; a band holding
    several swaps is bisected for each pair. The result is sorted and does not
    repeat the input stations.
    """
    zs = sorted(float(z) for z in zs)
    if len(zs) < 2:
        return ()
    tol = float(rel_tol) * max(abs(zs[-1] - zs[0]), 1.0)
    fracs = [_side_fracs_at(geo, z) for z in zs]
    found: set[float] = set()
    for (za, fa), (zb, fb) in zip(zip(zs, fracs), zip(zs[1:], fracs[1:])):
        if fa is None or fb is None:
            continue
        for side, n1, n2 in _swapped_pairs(fa, fb):
            z = _crossing_z(geo, side, n1, n2, za, zb, fa[side][n1] - fa[side][n2], tol)
            if z is not None and za < z < zb:
                found.add(round(z, 12))
    return tuple(sorted(found))


__all__ = ["key_crossing_stations"]
