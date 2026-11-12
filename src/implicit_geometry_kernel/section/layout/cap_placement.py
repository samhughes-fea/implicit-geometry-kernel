"""Declarative spar-cap placement rules for parametric SectionLayout."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

CapPlacementKind = Literal["none", "fixed_band", "per_web_pads", "continuous_box"]
FixedBandAnchor = Literal["pitching", "max_thickness"]


@dataclass(frozen=True)
class CapPlacementRule:
    """Cap placement policy (additive; SystemLayoutSpec remains catalogue truth)."""

    kind: CapPlacementKind
    anchor: FixedBandAnchor | None = None
    half_width: float | None = None
    chord_window: tuple[float, float] | None = None
    #: per-side chord windows ``((suction_lo, suction_hi), (pressure_lo, pressure_hi))``
    #: when the two caps differ; ``chord_window`` is then their union
    side_windows: tuple[tuple[float, float], tuple[float, float]] | None = None

    @staticmethod
    def none() -> CapPlacementRule:
        return CapPlacementRule(kind="none")

    @staticmethod
    def fixed_band(
        anchor: FixedBandAnchor = "max_thickness",
        half_width: float | None = None,
    ) -> CapPlacementRule:
        return CapPlacementRule(kind="fixed_band", anchor=anchor, half_width=half_width)

    @staticmethod
    def per_web_pads() -> CapPlacementRule:
        return CapPlacementRule(kind="per_web_pads")

    @staticmethod
    def continuous_box() -> CapPlacementRule:
        return CapPlacementRule(kind="continuous_box")

    @staticmethod
    def from_chord_windows(windows) -> CapPlacementRule:
        """Family-B band from published chord windows ``((lo, hi), ...)``.

        ``chord_window`` is the union interval in chord fraction from the LE.
        Assemble uses that interval, not a pitch-centred pad of ``half_width``.
        """
        pairs = [(float(a), float(b)) for a, b in (windows or ())]
        if not pairs:
            return CapPlacementRule.none()
        lo = min(min(p) for p in pairs)
        hi = max(max(p) for p in pairs)
        half = 0.5 * (hi - lo)
        if half <= 0.0:
            return CapPlacementRule.none()
        return CapPlacementRule(
            kind="fixed_band",
            anchor="pitching",
            half_width=half,
            chord_window=(lo, hi),
        )

    @staticmethod
    def from_side_windows(suction, pressure) -> CapPlacementRule:
        """Family-B caps with their own chord window on each side (fractions from the LE).

        Each ``(lo, hi)`` is where that side's cap is; ``chord_window`` is the union
        (for the wall-region graph's ordering), and the section builds each cap in its
        own window.
        """
        su = (float(min(suction)), float(max(suction)))
        pr = (float(min(pressure)), float(max(pressure)))
        if su[1] - su[0] <= 0.0 or pr[1] - pr[0] <= 0.0:
            raise ValueError(f"cap windows need lo < hi; got {suction} and {pressure}")
        rule = CapPlacementRule.from_chord_windows((su, pr))
        return CapPlacementRule(kind=rule.kind, anchor=rule.anchor, half_width=rule.half_width,
                                chord_window=rule.chord_window, side_windows=(su, pr))

    @staticmethod
    def from_structural_family(
        family: str,
        *,
        spar_cap_positioning: str | None = None,
    ) -> CapPlacementRule:
        fam = str(family).strip().upper()
        if fam == "A":
            return CapPlacementRule.none()
        if fam == "B":
            pos = str(spar_cap_positioning or "max_thickness").strip().lower()
            anchor: FixedBandAnchor = (
                "pitching" if pos in ("fixed", "pitching") else "max_thickness"
            )
            return CapPlacementRule.fixed_band(anchor=anchor)
        if fam == "C":
            return CapPlacementRule.per_web_pads()
        if fam == "D":
            return CapPlacementRule.continuous_box()
        raise ValueError(f"Unknown structural family {family!r}")


__all__ = [
    "CapPlacementKind",
    "CapPlacementRule",
    "FixedBandAnchor",
]