"""What a structure is made of: four kinds of contribution.

A thin-walled structure is a closed body (its outer mould line, OML) along a span,
whose wall may be partitioned and which may carry internal members. Builders turn
design parameters into contributions; :class:`~.assembler.StructureAssembler`
combines them. Every field is an :class:`~implicit_geometry_kernel.sdf.SDFNode`
evaluated at points ``(x, y, z)``; every length is in metres.

* :class:`WallContribution`: part of the wall. The base wall (no ``footprint``) covers
  the whole OML. A partition sets the wall's full thickness inside a footprint (one
  laminate, one mid-surface); a stacked layer instead lies beneath what is above it
  (a spar cap kept as its own layer under the skin).
* :class:`BranchContribution`: an internal member with its own mid-surface that meets
  the wall: shear webs, and stringers that stand off the wall with a free end. It may
  run over part of the span only.
* :class:`RegionContribution`: a volume with no mid-surface (core, bond).
* :class:`TransverseContribution`: a member lying in the section plane at one station
  (rib, bulkhead, frame).

A new kind of member that fits one of these needs only a new builder. Every
contribution may carry a ``source`` (whatever object made it, such as a section
component); the assembler never reads it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from implicit_geometry_kernel.sdf.nodes.base import FeatureSpec, SDFNode
from implicit_geometry_kernel.structure.laws import SpanLaw


@dataclass(frozen=True)
class WallContribution:
    """Part of the wall.

    ``footprint`` is negative where the contribution applies; ``None`` marks the base
    wall. For the base wall and a partition, ``thickness`` is the full wall thickness
    where it applies; where partitions overlap, the higher priority sets it. For a
    ``stacked`` layer, ``thickness`` is the layer's own thickness, added beneath the
    wall above it (layers stack in contribution order).
    """

    spec: FeatureSpec
    thickness: SpanLaw
    footprint: SDFNode | None = None
    stacked: bool = False
    source: Any = field(default=None, compare=False, repr=False)

    def __post_init__(self) -> None:
        if self.stacked and self.footprint is None:
            raise ValueError(f"{self.spec.name}: a stacked layer needs a footprint")

    @property
    def is_base(self) -> bool:
        return self.footprint is None


@dataclass(frozen=True)
class BranchContribution:
    """An internal member: mid-surface ``mid`` and material ``region``.

    ``clip`` (negative where allowed) restricts where the mid-surface is valid, for
    example to one side of the chord line. ``depth`` bounds how far the member stands
    off the wall mid-surface; a member with a depth has a ``free_end``. ``extent`` is
    the span interval ``(z_start, z_end)`` it occupies (``None``: the whole span); its
    ends are topology events, which meshes should treat as forced stations.
    """

    spec: FeatureSpec
    mid: SDFNode
    region: SDFNode
    thickness: SpanLaw
    clip: SDFNode | None = None
    depth: SpanLaw | None = None
    free_end: bool = False
    extent: tuple[float, float] | None = None
    source: Any = field(default=None, compare=False, repr=False)

    def __post_init__(self) -> None:
        if self.extent is not None and not self.extent[1] >= self.extent[0]:
            raise ValueError(f"{self.spec.name}: extent must be (z_start, z_end) with z_end >= z_start")
        if self.free_end and self.depth is None:
            raise ValueError(f"{self.spec.name}: a free end needs a depth")

    def present(self, z: float, tol: float = 1e-9) -> bool:
        """True when the member exists at station ``z``."""
        return self.extent is None or self.extent[0] - tol <= z <= self.extent[1] + tol


@dataclass(frozen=True)
class RegionContribution:
    """A volume with no mid-surface. ``tag_region`` (default ``region``) selects what
    downstream tools should label as belonging to it."""

    spec: FeatureSpec
    region: SDFNode
    tag_region: SDFNode | None = None
    source: Any = field(default=None, compare=False, repr=False)

    @property
    def tag_field(self) -> SDFNode:
        return self.region if self.tag_region is None else self.tag_region


@dataclass(frozen=True)
class TransverseContribution:
    """A member lying in the section plane at station ``z`` (rib, bulkhead, frame).

    ``region`` is its in-plane shape (negative inside); ``thickness`` its spanwise
    thickness. It is not part of any section cut; it sits between them.
    """

    spec: FeatureSpec
    z: float
    region: SDFNode
    thickness: float
    source: Any = field(default=None, compare=False, repr=False)

    def __post_init__(self) -> None:
        if not self.thickness > 0.0:
            raise ValueError(f"{self.spec.name}: thickness must be positive")


Contribution = WallContribution | BranchContribution | RegionContribution | TransverseContribution

__all__ = [
    "BranchContribution", "Contribution", "RegionContribution", "TransverseContribution",
    "WallContribution",
]
