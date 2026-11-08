"""Named skin-loop partitions that need a constitutive laminate OFFSET.

The kernel does not emit ABD. Downstream constitutives apply eccentricity only
inside these regions. Mesh xy stays the skin mid everywhere.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LaminateOffsetRegion:
    """Cap (or later doubler) footprint on the conformal skin mid.

    ``eccentricity_m`` is the equal-modulus stack mid in the inward shell
    normal ``+n`` from ``n = 0`` at the skin mid (``t_cap / 2`` for a
    ``(cap, skin)`` stack). Stiffness-weighted NA is the consumer's job.
    """

    name: str
    requires_eccentricity: bool
    eccentricity_m: float
    layers: tuple[str, ...]
    layer_thicknesses_m: tuple[float, ...]
    skin_side: str | None = None
    reference_surface: str = "skin_mid"

    @property
    def id(self) -> str:
        return self.name


__all__ = ["LaminateOffsetRegion"]
