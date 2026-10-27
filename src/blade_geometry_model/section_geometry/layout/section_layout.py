"""Parametric SectionLayout facade over the SystemType catalogue."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from blade_geometry_model.section_geometry.layout.cap_placement import (
    CapPlacementRule,
)
from blade_geometry_model.section_geometry.layout.system_layout import (
    SystemLayoutSpec,
    WebOrientation,
    resolve_system_type,
)

LayoutFlag = Literal["has_spar_caps", "zero_web_catalogue"]


@dataclass(frozen=True)
class SectionLayout:
    """Parametric layout: web fractions + orientation + cap placement.

    Absolute web stations (``web_stations_m``) are optional overrides; catalogue
    conversion stores chord fractions. Catalogue keys are **presets** —
    :func:`sections.from_layout.multicell_from_section_layout` assembles without
    requiring ``catalogue_key``. ``to_system_layout_spec`` remains best-effort
    round-trip when a key is set.
    """

    web_chord_fracs: tuple[float, ...]
    orientation: WebOrientation
    cap_placement: CapPlacementRule
    catalogue_key: str | None = None
    web_stations_m: tuple[float, ...] | None = None
    # Optional non-parallel axes: ((x_top, y_top), (x_bot, y_bot)) per web [m].
    # When set, chord fractions / stations are ignored for assemble geometry.
    web_axes_m: (
        tuple[tuple[tuple[float, float], tuple[float, float]], ...] | None
    ) = None
    structural_family: str | None = None
    flags: tuple[str, ...] = field(default_factory=tuple)

    @property
    def n_webs(self) -> int:
        if self.web_axes_m is not None:
            return len(self.web_axes_m)
        if self.web_stations_m is not None:
            return len(self.web_stations_m)
        return len(self.web_chord_fracs)

    @property
    def key(self) -> str | None:
        """Catalogue preset id when this layout came from (or re-attached to) the registry."""
        return self.catalogue_key

    @property
    def has_spar_caps(self) -> bool:
        """True when cap placement is not ``none`` (or the catalogue flag is set)."""
        if "has_spar_caps" in self.flags:
            return True
        return self.cap_placement.kind != "none"

    @property
    def web_orientation(self) -> WebOrientation:
        """Alias for :attr:`orientation` (SystemLayoutSpec vocabulary)."""
        return self.orientation

    @staticmethod
    def from_catalogue_key(key: str) -> SectionLayout:
        return SectionLayout.from_system_layout_spec(resolve_system_type(key))

    @staticmethod
    def from_system_layout_spec(spec: SystemLayoutSpec) -> SectionLayout:
        """Lift a catalogue preset into the parametric :class:`SectionLayout` SoT."""
        cap = CapPlacementRule.from_structural_family(
            spec.structural_family,
            spar_cap_positioning=spec.spar_cap_positioning,
        )
        flags: list[str] = []
        if spec.has_spar_caps:
            flags.append("has_spar_caps")
        if spec.geometry_mode == "zero_web_catalogue":
            flags.append("zero_web_catalogue")
        return SectionLayout(
            web_chord_fracs=tuple(float(x) for x in spec.web_chord_fracs),
            orientation=spec.web_orientation,
            cap_placement=cap,
            catalogue_key=spec.key,
            web_stations_m=None,
            web_axes_m=None,
            structural_family=str(spec.structural_family),
            flags=tuple(flags),
        )

    def to_system_layout_spec(self) -> SystemLayoutSpec:
        """Best-effort catalogue round-trip when ``catalogue_key`` is set."""
        if not self.catalogue_key:
            raise ValueError(
                "to_system_layout_spec requires catalogue_key; "
                "SystemLayoutSpec remains the registry source of truth."
            )
        return resolve_system_type(self.catalogue_key)

    def with_web_stations_m(self, stations_m: tuple[float, ...]) -> SectionLayout:
        return SectionLayout(
            web_chord_fracs=self.web_chord_fracs,
            orientation=self.orientation,
            cap_placement=self.cap_placement,
            catalogue_key=self.catalogue_key,
            web_stations_m=tuple(float(x) for x in stations_m),
            web_axes_m=self.web_axes_m,
            structural_family=self.structural_family,
            flags=self.flags,
        )

    def to_jsonable(self) -> dict[str, Any]:
        return {
            "web_chord_fracs": list(self.web_chord_fracs),
            "orientation": self.orientation,
            "cap_placement": {
                "kind": self.cap_placement.kind,
                "anchor": self.cap_placement.anchor,
                "half_width": self.cap_placement.half_width,
            },
            "catalogue_key": self.catalogue_key,
            "web_stations_m": (
                None if self.web_stations_m is None else list(self.web_stations_m)
            ),
            "web_axes_m": (
                None
                if self.web_axes_m is None
                else [
                    [[float(a), float(b)] for a, b in axis] for axis in self.web_axes_m
                ]
            ),
            "structural_family": self.structural_family,
            "flags": list(self.flags),
        }


__all__ = [
    "LayoutFlag",
    "SectionLayout",
]