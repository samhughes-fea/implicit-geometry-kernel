"""Hard-coded web positions for catalogue LC / LE / PA anchoring.

LC and LE store authored positions in their native units (chord fraction; metres
from the leading edge). PA stores blade-x from the pitch axis: primary webs are
uniform interiors of the root ``[x_LE, x_TE]`` band; extra webs continue that
same lattice (first extra at root TE). Those webs run parallel to the pitch axis.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

_MAX_PRIMARY = 3


@dataclass(frozen=True)
class PaRegions:
    """Planform blade-x bands for PA webs (metres from the pitch axis)."""

    x_le_root: float
    x_te_root: float
    x_te_max_chord: float
    z_max_chord: float

    @property
    def primary(self) -> tuple[float, float]:
        lo, hi = (self.x_le_root, self.x_te_root)
        return (lo, hi) if lo <= hi else (hi, lo)

    @property
    def secondary(self) -> tuple[float, float]:
        lo, hi = (self.x_te_root, self.x_te_max_chord)
        return (lo, hi) if lo <= hi else (hi, lo)

    @property
    def root_half_width(self) -> float:
        return float(min(abs(self.x_le_root), abs(self.x_te_root)))

    @property
    def secondary_width(self) -> float:
        return float(self.x_te_max_chord - self.x_te_root)


def chord_fraction_from_pitch(x_m: float, chord: float, pivot: float) -> float:
    """Chord fraction of a point ``x_m`` metres from the pitch axis, toward the TE.

    ``pivot`` is the pitch-axis fraction of this station's own chord. The same
    conversion places webs, cap edges and stringers that were measured from the
    pitch axis on the chord law.
    """
    c = float(chord)
    if c <= 0.0:
        raise ValueError("chord must be positive")
    return float(pivot) + float(x_m) / c


def blade_x_le_te(chord: float, pitch: float) -> tuple[float, float]:
    """Blade-x of the chord LE and TE at a station (pitch axis at ``x = 0``)."""
    c, p = float(chord), float(pitch)
    return (0.0 - p) * c, (1.0 - p) * c


def pa_regions(
    stations: Sequence,
    pitch_at: Callable[[float], float],
) -> PaRegions:
    """Primary = root LE/TE; secondary = root TE to max-chord TE."""
    if not stations:
        raise ValueError("pa_regions needs at least one station")
    root = stations[0]
    p0 = float(pitch_at(float(root.z)))
    x_le0, x_te0 = blade_x_le_te(float(root.chord), p0)
    mc = max(stations, key=lambda s: float(s.chord))
    pmc = float(pitch_at(float(mc.z)))
    _le_mc, x_te_mc = blade_x_le_te(float(mc.chord), pmc)
    return PaRegions(
        x_le_root=float(x_le0),
        x_te_root=float(x_te0),
        x_te_max_chord=float(x_te_mc),
        z_max_chord=float(mc.z),
    )


def pa_web_x_from_pitch(
    stations: Sequence,
    pitch_at: Callable[[float], float],
    n_webs: int,
) -> tuple[float, ...]:
    """Uniform PA blade-x list of length ``n_webs``.

    Up to three primary webs sit at interiors of the root ``[x_LE, x_TE]``
    band (spacing ``Δx = (x_TE − x_LE) / (n_primary + 1)``). Extra webs
    continue that lattice: slot ``n_primary + k`` at ``x_LE + (n_primary + k)·Δx``
    (first extra at root TE).
    """
    n = int(n_webs)
    if n <= 0:
        return ()
    regions = pa_regions(stations, pitch_at)
    x_min, x_max = float(regions.x_le_root), float(regions.x_te_root)
    n_primary = min(n, _MAX_PRIMARY)
    n_secondary = n - n_primary
    dx = (x_max - x_min) / (n_primary + 1)
    return tuple(x_min + i * dx for i in range(1, n_primary + n_secondary + 1))


def lc_web_fractions(fracs: Sequence[float]) -> tuple[float, ...]:
    """Authored LC chord fractions (hard-coded; remapped into local chord)."""
    return tuple(float(f) for f in fracs)


def le_web_offsets_m(fracs: Sequence[float], root_chord: float) -> tuple[float, ...]:
    """Authored LE offsets: catalogue fractions times the root chord (metres)."""
    c0 = float(root_chord)
    return tuple(float(f) * c0 for f in fracs)


__all__ = [
    "PaRegions",
    "blade_x_le_te",
    "chord_fraction_from_pitch",
    "lc_web_fractions",
    "le_web_offsets_m",
    "pa_regions",
    "pa_web_x_from_pitch",
]
