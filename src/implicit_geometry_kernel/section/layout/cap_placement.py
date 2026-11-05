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