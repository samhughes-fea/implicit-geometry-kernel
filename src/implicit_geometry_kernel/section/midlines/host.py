"""Public section host Protocol for midline / mesh export (C++-portable seam)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Protocol, cast, runtime_checkable


@runtime_checkable
class ShellComponentLike(Protocol):
    """Chord-frame shell solid that can emit a midline polyline.

    Narrower Protocols
    (:class:`SkinShellComponentLike`, :class:`WebShellComponentLike`,
    :class:`CapShellComponentLike`) document the ABI shape of thickness /
    alignment / cap accessors. Strip **kind** is chosen by recipe +
    :func:`~implicit_geometry_kernel.section.assembly.members.wall_type_for_label`,
    not by ``isinstance`` on these Protocols.
    """

    def midline_polyline(self, n: int = ...) -> Any:
        ...


@runtime_checkable
class SkinShellComponentLike(Protocol):
    """Outer-skin solid ABI: closed midline + wall thickness.

    Kind discrimination stays on ``wall_type_for_label`` / recipes.
    """

    thickness: float

    def midline_polyline(self) -> Any:
        ...


@runtime_checkable
class WebShellComponentLike(Protocol):
    """Shear-web solid ABI: open midline + thickness + alignment.

    Kind discrimination stays on ``wall_type_for_label`` / recipes.
    """

    thickness: float
    alignment: str

    def midline_polyline(self, n: int = ...) -> Any:
        ...


@runtime_checkable
class CapShellComponentLike(Protocol):
    """Spar-cap solid ABI: open midline + ``cap_height`` + ``surface``.

    Solids emit via :func:`~implicit_geometry_kernel.section.midlines.export.build_solid_cap_strips`
    reading recipe attrs; kind is not decided by ``isinstance`` on this Protocol.
    """

    cap_height: float
    surface: str

    def midline_polyline(self, n: int = ...) -> Any:
        ...


@runtime_checkable
class SectionMidlineHost(Protocol):
    """Section objects that can emit FE-agnostic midline strips.

    Implementations expose chord-frame component maps via *public* properties
    (not private ``_`` attributes) so mesh / export seams stay ABI-shaped.
    This is the DIP seam; Skin/Web/Cap Protocols below are shape docs only.
    """
    @property
    def airfoil(self) -> Any:
        """Airfoil SDF / outline carrier used for LE/TE and chord."""
        ...

    @property
    def twist_rad(self) -> float:
        """Section twist angle in radians (S → B frame)."""
        ...

    @property
    def components_unrotated(self) -> Mapping[str, Any]:
        """Chord-frame skin + web (+ core) components keyed by label."""
        ...

    @property
    def spar_cap_components_unrotated(self) -> Mapping[str, Any]:
        """Chord-frame spar-cap components keyed by label (may be empty)."""
        ...

    @property
    def components(self) -> Mapping[str, Any]:
        """Active labelled component map. May be twisted; emit prefers the chord-frame map."""
        ...

    @property
    def labels(self) -> Any:
        """Component labels in map order."""
        ...

    def __getitem__(self, key: str) -> Any:
        """Component stored under ``key``."""
        ...


def resolve_section_host(section: Any) -> SectionMidlineHost:
    """Unwrap wrappers that expose ``underlying_section`` when present."""
    host = section.underlying_section if hasattr(section, "underlying_section") else section
    return cast(SectionMidlineHost, host)


__all__ = [
    "ShellComponentLike",
    "SkinShellComponentLike",
    "WebShellComponentLike",
    "CapShellComponentLike",
    "SectionMidlineHost",
    "resolve_section_host",
]
