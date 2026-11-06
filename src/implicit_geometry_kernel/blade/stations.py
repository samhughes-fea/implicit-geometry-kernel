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
from typing import Tuple

from implicit_geometry_kernel.blade.airfoil import Airfoil


# ---------------------------------------------------------------------------
# Blade definition
# ---------------------------------------------------------------------------
@dataclass
class Station:
    """A spanwise station: position ``z``, chord, twist (degrees, nose up about +z) and airfoil.

    ``pitch_axis`` is the chord fraction that lies on the blade axis (x = 0), and
    ``offset`` an (x, y) shift of that point (prebend, sweep).
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


def __getattr__(name: str):
    if name == "SpanWeb":
        import warnings

        warnings.warn(
            "SpanWeb is deprecated; use ShearWeb. Removed in 0.6.0.",
            DeprecationWarning,
            stacklevel=2,
        )
        return ShearWeb
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
