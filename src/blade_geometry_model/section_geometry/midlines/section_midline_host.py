"""Public section host Protocol for midline / mesh export (C++-portable seam)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class ShellComponentLike(Protocol):
    """Chord-frame shell solid that can emit a midline polyline.

    Narrower runtime-checkable Protocols
    (:class:`SkinShellComponentLike`, :class:`WebShellComponentLike`,
    :class:`CapShellComponentLike`) add the thickness / alignment / cap
    accessors that :mod:`shell_midline_export` uses for strip kind.
    """

    def midline_polyline(self, n: int = ...) -> Any:
        ...


@runtime_checkable
class SkinShellComponentLike(Protocol):
    """Outer-skin solid: closed midline + wall thickness."""

    thickness: float

    def midline_polyline(self) -> Any:
        ...


@runtime_checkable
class WebShellComponentLike(Protocol):
    """Shear-web solid: open midline + thickness + alignment."""

    thickness: float
    alignment: str

    def midline_polyline(self, n: int = ...) -> Any:
        ...


@runtime_checkable
class CapShellComponentLike(Protocol):
    """Spar-cap solid: open midline + cap height + surface side."""

    cap_height: float
    surface: str

    def midline_polyline(self, n: int = ...) -> Any:
        ...


@runtime_checkable
class SectionMidlineHost(Protocol):
    """Section objects that can emit FE-agnostic midline strips.

    Implementations expose chord-frame component maps via *public* properties
    (not private ``_`` attributes) so mesh / export seams stay ABI-shaped.
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


def resolve_section_host(section: Any) -> Any:
    """Unwrap wrappers that expose ``underlying_section`` when present."""
    if hasattr(section, "underlying_section"):
        return section.underlying_section
    return section


__all__ = [
    "ShellComponentLike",
    "SkinShellComponentLike",
    "WebShellComponentLike",
    "CapShellComponentLike",
    "SectionMidlineHost",
    "resolve_section_host",
]
