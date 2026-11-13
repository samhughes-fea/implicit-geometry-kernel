"""Blade definition: spanwise stations, skin regions and shear webs.

These types describe the shape and the structural layout of a blade or wing, and
nothing about materials: ``Station`` (span position, chord, twist, pitch axis and
offset), ``SkinRegion`` (a chordwise band of skin on both surfaces) and ``ShearWeb``
(a web at a chord fraction). A downstream shell mesher turns them into shells.

Axes: z runs along the span, x along the chord (leading edge to trailing edge) and
y through the thickness (suction side positive). The pitch axis sits at x = 0.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from implicit_geometry_kernel.blade.airfoil import Airfoil


# ---------------------------------------------------------------------------
# Blade definition
# ---------------------------------------------------------------------------
@dataclass
class Station:
    """A spanwise station: position ``z``, chord, twist (degrees, nose up about +z) and airfoil.

    ``pitch_axis`` is stored for CSV / windIO round-trip. Loft stacking and the
    S-frame origin use required ``BladeGeometry.pivot_fraction``, not this default.
    ``offset`` is an (x, y) shift of that piercing point (prebend, sweep).
    """

    z: float
    chord: float
    airfoil: Airfoil
    twist: float = 0.0
    pitch_axis: float = 0.25
    offset: Tuple[float, float] = (0.0, 0.0)

    def __post_init__(self):
        if self.chord <= 0.0:
            raise ValueError("chord must be positive")


@dataclass(frozen=True)
class TransverseMember:
    """A rib, bulkhead or frame lying in the section plane at one span station.

    ``region`` is the cavity the face fills (cell interiors of the cut, or the
    skin loop when there are no cells). Not a wall type and not a catalogue letter.
    """

    name: str
    z: float
    thickness: float
    region: str = "cavity"

    def __post_init__(self):
        if not str(self.name).strip():
            raise ValueError("name is required")
        if self.thickness <= 0.0:
            raise ValueError("thickness must be positive")
        if not str(self.region).strip():
            raise ValueError("region is required")


@dataclass(frozen=True)
class TipClosure:
    """The face that closes the skin shell at the blade tip (z/L = 1).

    Every ``BladeGeometry`` has one. It is part of the skin, not a rib: its faces
    carry the skin partition (``partition``) and so the skin's material, and it is
    not counted in a layout's ``nR``. It fills the skin mid-surface loop at the tip
    station, split at any web that reaches the tip. ``thickness`` of ``None`` takes
    the skin thickness law at the tip.
    """

    name: str = "tip_closure"
    thickness: Optional[float] = None
    partition: str = "skin"

    def __post_init__(self):
        if not str(self.name).strip():
            raise ValueError("name is required")
        if self.thickness is not None and self.thickness <= 0.0:
            raise ValueError("thickness must be positive")
        if not str(self.partition).strip():
            raise ValueError("partition is required")


@dataclass(frozen=True)
class SkinRegion:
    """Skin between chord fractions ``start`` and ``end`` on both surfaces, with ``n`` elements."""

    name: str
    start: float
    end: float
    n: int = 4

    def __post_init__(self):
        if not 0.0 <= self.start < self.end <= 1.0:
            raise ValueError(f"region '{self.name}' needs 0 <= start < end <= 1")
        if self.n < 1:
            raise ValueError("a region needs at least one element")


@dataclass(frozen=True)
class ShearWeb:
    """A shear web joining the surfaces at chord fraction ``position``, meshed with ``n`` shell elements."""

    name: str
    position: float
    n: int = 4


@dataclass(frozen=True)
class Stringer:
    """A free-end or through branch standing off one skin, optionally over part of the span.

    ``extent`` is the inclusive ``(z_start, z_end)`` interval; ``None`` means the
    whole blade. Ends are topology events and ``forced_stations()`` mesh hints.

    ``position`` is a chord fraction from the leading edge. ``position_m`` is
    metres from the pitch axis toward the trailing edge; when it is set, the
    cut converts it with that station's chord and pivot and ``position`` is
    only a placeholder.
    """

    name: str
    position: float
    thickness: float
    depth: float
    extent: Optional[Tuple[float, float]] = None
    side: str = "pressure"
    free_end: bool = True
    position_m: Optional[float] = None

    def __post_init__(self):
        if not str(self.name).strip():
            raise ValueError("name is required")
        if self.position_m is None and not 0.0 < self.position < 1.0:
            raise ValueError("position must be a chord fraction inside (0, 1)")
        if self.position_m is not None:
            object.__setattr__(self, "position_m", float(self.position_m))
        if self.thickness <= 0.0:
            raise ValueError("thickness must be positive")
        if self.depth <= 0.0:
            raise ValueError("depth must be positive")
        if self.side not in ("pressure", "suction"):
            raise ValueError("side must be 'pressure' or 'suction'")
        if self.extent is not None:
            z0, z1 = float(self.extent[0]), float(self.extent[1])
            if z0 >= z1:
                raise ValueError("extent needs z_start < z_end")
            object.__setattr__(self, "extent", (z0, z1))


@dataclass(frozen=True)
class PlyDrop:
    """Geometry-only ply-drop ramp: extra thickness inside a (z, s) window.

    ``s`` is chord fraction from the leading edge. Thickness ramps linearly
    over ``taper`` in both s and z. Not a layup library.
    """

    name: str
    s0: float
    s1: float
    z0: float
    z1: float
    thickness: float
    taper: float = 0.05

    def __post_init__(self):
        if not str(self.name).strip():
            raise ValueError("name is required")
        if not 0.0 <= self.s0 < self.s1 <= 1.0:
            raise ValueError("ply-drop needs 0 <= s0 < s1 <= 1")
        if self.z0 >= self.z1:
            raise ValueError("ply-drop needs z0 < z1")
        if self.thickness < 0.0 or self.taper < 0.0:
            raise ValueError("thickness and taper must be >= 0")

    def extra(self, z: float, s: float) -> float:
        wz = _ramp(float(z), self.z0, self.z1, self.taper)
        ws = _ramp(float(s), self.s0, self.s1, self.taper)
        return float(self.thickness) * min(wz, ws)


def _ramp(x: float, lo: float, hi: float, taper: float) -> float:
    if x < lo - taper or x > hi + taper:
        return 0.0
    if x < lo:
        return 0.0 if taper <= 0.0 else (x - (lo - taper)) / taper
    if x > hi:
        return 0.0 if taper <= 0.0 else ((hi + taper) - x) / taper
    return 1.0
