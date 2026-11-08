"""Hard-coded web positions for catalogue LC / LE / PA anchoring.

LC and LE store authored positions in their native units (chord fraction; metres
from the leading edge). PA stores blade-x from the pitch axis and, uniquely,
lives in two planform regions because those webs run parallel to the pitch axis.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

# Authored PA slots. ``n_webs`` takes a prefix (``n == 1`` is the PA itself).
# Primary ``t`` is a signed fraction of the root half-width toward PA.
# Secondary ``t`` is a fraction along [x_TE(root), x_TE(max_chord)].
_PA_SLOTS: tuple[tuple[str, float], ...] = (
    ("primary", -0.35),
    ("primary", 0.35),
    ("primary", 0.0),
    ("secondary", 0.40),
    ("secondary", 0.70),
)


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


def _slots_for(n_webs: int) -> tuple[tuple[str, float], ...]:
    n = int(n_webs)
    if n <= 0:
        return ()
    if n == 1:
        return (("primary", 0.0),)
    if n <= len(_PA_SLOTS):
        return _PA_SLOTS[:n]
    extra: list[tuple[str, float]] = []
    k = n - len(_PA_SLOTS)
    for i in range(k):
        extra.append(("secondary", (i + 1) / (k + 1)))
    return _PA_SLOTS + tuple(extra)


def _in_interval(x: float, bounds: tuple[float, float], *, tol: float = 1e-9) -> bool:
    lo, hi = bounds
    return float(lo) - tol <= float(x) <= float(hi) + tol


def pa_web_x_from_pitch(
    stations: Sequence,
    pitch_at: Callable[[float], float],
    n_webs: int,
) -> tuple[float, ...]:
    """Hard-coded PA blade-x list of length ``n_webs``, validated against regions."""
    n = int(n_webs)
    if n <= 0:
        return ()
    regions = pa_regions(stations, pitch_at)
    half = max(regions.root_half_width, 1e-12)
    shelf = regions.secondary_width
    xs: list[float] = []
    for kind, t in _slots_for(n):
        if kind == "primary":
            x = float(t) * half
            if not _in_interval(x, regions.primary):
                raise ValueError(
                    f"PA primary web x={x} is outside root bounds {regions.primary}"
                )
            xs.append(x)
            continue
        if shelf <= 1e-9:
            x = float(t) * half * 0.5
            if not _in_interval(x, regions.primary):
                raise ValueError(
                    f"PA extra web x={x} is outside root bounds {regions.primary} "
                    "(no secondary shelf)"
                )
            xs.append(x)
            continue
        x = float(regions.x_te_root) + float(t) * shelf
        if not _in_interval(x, regions.secondary):
            raise ValueError(
                f"PA secondary web x={x} is outside shelf {regions.secondary}"
            )
        xs.append(x)
    return tuple(xs)


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
    "lc_web_fractions",
    "le_web_offsets_m",
    "pa_regions",
    "pa_web_x_from_pitch",
]
