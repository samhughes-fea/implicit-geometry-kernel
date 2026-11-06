"""Combine contributions into one structure.

The wall is described by the OML and a thickness field ``T(x, y, z)``: the base wall
thickness, replaced by a partition's thickness inside its footprint, plus any stacked
layers beneath. From these:

* the wall reference (mid) surface is ``phi_OML + T / 2``, a distance-like field whose
  zero set lies half the local wall thickness inside the OML;
* the cavity is ``phi_OML + T`` (negative in the hollow inside the wall);
* the base wall and partitions own the band from the OML down to the partitions'
  thickness (inside its footprint for a partition); each stacked layer owns its own
  band beneath that.

Branches are valid inside the wall reference surface, within ``depth`` of it when they
have one, and inside their ``clip``. Every material region joins a priority union
(higher ``FeatureSpec.priority`` owns overlaps), hard-clipped to the OML.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from implicit_geometry_kernel.sdf.field import EvalCache, FeatureSpec, FeatureTable, SDFNode
from implicit_geometry_kernel.sdf.nodes import FieldOffsetNode, PriorityUnionNode, TagNode, TrimNode
from implicit_geometry_kernel.sdf.span_nodes import ComplementNode, LawOffsetNode, MaxNode, SpanIntervalNode
from implicit_geometry_kernel.structure.builders import BuildContext, BuilderRegistry
from implicit_geometry_kernel.structure.contributions import (
    BranchContribution,
    Contribution,
    RegionContribution,
    TransverseContribution,
    WallContribution,
)


class StructureError(ValueError):
    """The contributions do not describe a valid structure."""


class WallThickness:
    """``T(points)``: the wall thickness at each point.

    Start from the base wall; inside each partition's footprint the partition's
    thickness replaces it (partitions applied in rising priority, so the highest
    wins); then each stacked layer adds its own thickness inside its footprint.
    """

    def __init__(self, walls: tuple[WallContribution, ...]) -> None:
        bases = [w for w in walls if w.is_base]
        if len(bases) != 1:
            raise StructureError(f"the wall needs exactly one base contribution, got {len(bases)}")
        self.base = bases[0]
        self.partitions = tuple(sorted((w for w in walls if not w.is_base and not w.stacked),
                                       key=lambda w: w.spec.priority))
        self.layers = tuple(w for w in walls if w.stacked)

    def __call__(self, points: np.ndarray) -> np.ndarray:
        return self.up_to(len(self.layers))(points)

    def up_to(self, n_layers: int) -> _PartialThickness:
        """Thickness of the wall above stacked layer ``n_layers`` (0: partitions only)."""
        return _PartialThickness(self, n_layers)

    def _evaluate(self, points: np.ndarray, n_layers: int) -> np.ndarray:
        pts = np.asarray(points, dtype=np.float64)
        z = pts[:, 2]
        t = np.asarray(self.base.thickness(z), dtype=np.float64).copy()
        cache = EvalCache(pts)
        for p in self.partitions:
            inside = p.footprint.evaluate(pts, cache).phi < 0.0
            if np.any(inside):
                t = np.where(inside, p.thickness(z), t)
        for layer in self.layers[:n_layers]:
            inside = layer.footprint.evaluate(pts, cache).phi < 0.0
            if np.any(inside):
                t = t + np.where(inside, layer.thickness(z), 0.0)
        return t


class _PartialThickness:
    """Picklable ``points -> T`` counting only the first ``n_layers`` stacked layers."""

    def __init__(self, wall: WallThickness, n_layers: int) -> None:
        self.wall, self.n_layers = wall, int(n_layers)

    def __call__(self, points: np.ndarray) -> np.ndarray:
        return self.wall._evaluate(points, self.n_layers)


@dataclass(frozen=True)
class BranchField:
    """A branch's mid-surface restricted to ``mask < 0``."""

    contribution: BranchContribution
    mask: SDFNode

    @property
    def spec(self) -> FeatureSpec:
        return self.contribution.spec


@dataclass(frozen=True)
class Structure:
    """The assembled structure: every field a downstream tool needs.

    ``material`` is negative inside any member and carries the owning member's id as
    ``meta``. ``wall_mid`` and ``cavity`` follow the module notes above.
    """

    table: FeatureTable
    span: tuple[float, float]
    oml: SDFNode
    thickness: WallThickness
    wall_mid: SDFNode
    cavity: SDFNode
    wall_meta: SDFNode
    material: SDFNode
    walls: tuple[WallContribution, ...]
    branches: tuple[BranchField, ...]
    regions: tuple[RegionContribution, ...]
    transverse: tuple[TransverseContribution, ...]

    def _interior_ends(self, i: int) -> tuple[float, ...]:
        z0, z1 = self.span
        return tuple(sorted({
            b.contribution.extent[i] for b in self.branches
            if b.contribution.extent is not None and z0 < b.contribution.extent[i] < z1
        }))

    @property
    def member_starts(self) -> tuple[float, ...]:
        """Stations inside the span where a branch starts."""
        return self._interior_ends(0)

    @property
    def member_ends(self) -> tuple[float, ...]:
        """Stations inside the span where a branch ends."""
        return self._interior_ends(1)

    @property
    def events(self) -> tuple[float, ...]:
        """Stations where the section topology changes: branch starts and ends."""
        return tuple(sorted({*self.member_starts, *self.member_ends}))

    @property
    def transverse_stations(self) -> tuple[float, ...]:
        """Stations of transverse members, in span order."""
        return tuple(t.z for t in self.transverse)

    def forced_stations(self) -> tuple[float, ...]:
        """Every station a spanwise mesh must contain: events and transverse members."""
        return tuple(sorted({*self.events, *self.transverse_stations}))

    def members_at(self, z: float) -> tuple[str, ...]:
        """Names of the wall, branch and region members present in the cut at ``z``."""
        names = [w.spec.name for w in self.walls]
        names += [b.spec.name for b in self.branches if b.contribution.present(z)]
        names += [r.spec.name for r in self.regions]
        return tuple(names)


class StructureAssembler:
    """Runs the builders of a registry and assembles the result."""

    def __init__(self, registry: BuilderRegistry) -> None:
        self.registry = registry

    def assemble(self, ctx: BuildContext) -> Structure:
        return assemble(ctx, self.registry.run(ctx))


def _check_ids(table: FeatureTable, contributions: list[Contribution]) -> None:
    seen: set[int] = set()
    for c in contributions:
        if c.spec.id in seen:
            raise StructureError(f"feature {c.spec.name!r} contributes twice")
        seen.add(c.spec.id)
        if c.spec.id >= len(table) or table[c.spec.id] != c.spec:
            raise StructureError(f"feature {c.spec.name!r} is not in the build table")


def assemble(ctx: BuildContext, contributions: list[Contribution]) -> Structure:
    """Combine ``contributions`` (all made with ``ctx.table``) into a :class:`Structure`."""
    _check_ids(ctx.table, contributions)
    oml = ctx.oml
    walls = tuple(c for c in contributions if isinstance(c, WallContribution))
    if not walls:
        raise StructureError("the structure needs at least one wall contribution")
    thickness = WallThickness(walls)
    wall_mid = FieldOffsetNode(oml, thickness, 0.5)
    cavity = FieldOffsetNode(oml, thickness, 1.0)
    # Base wall and partitions fill from the OML down to the partitions-only depth;
    # stacked layer k fills from depth T_k to T_(k+1) inside its footprint.
    envelope = ComplementNode(FieldOffsetNode(oml, thickness.up_to(0), 1.0))
    wall_region: dict[int, SDFNode] = {}
    for w in walls:
        if w.stacked:
            k = thickness.layers.index(w)
            above = FieldOffsetNode(oml, thickness.up_to(k), 1.0)
            below = ComplementNode(FieldOffsetNode(oml, thickness.up_to(k + 1), 1.0))
            r = MaxNode(above, below, w.footprint, meta_from=0)
        elif w.is_base:
            # A partition replaces the base wall inside its footprint.
            outside = [ComplementNode(p.footprint) for p in thickness.partitions]
            r = MaxNode(envelope, *outside, meta_from=0) if outside else envelope
        else:
            r = MaxNode(envelope, w.footprint, meta_from=0)
        wall_region[w.spec.id] = TagNode(TrimNode(r, oml), w.spec.id)
    wall_meta = PriorityUnionNode([wall_region[w.spec.id] for w in walls],
                                  [w.spec.priority for w in walls])

    branches: list[BranchField] = []
    material: list[SDFNode] = []
    priority: list[int] = []
    for c in contributions:
        if isinstance(c, WallContribution):
            region = wall_region[c.spec.id]
        elif isinstance(c, BranchContribution):
            masks: list[SDFNode] = [wall_mid]
            if c.depth is not None:
                masks.append(ComplementNode(LawOffsetNode(wall_mid, c.depth, scale=-1.0)))
            if c.clip is not None:
                masks.append(c.clip)
            if c.extent is not None:
                masks = [SpanIntervalNode(m, *c.extent) for m in masks[:1]] + masks[1:]
            branches.append(BranchField(c, masks[0] if len(masks) == 1 else MaxNode(*masks, meta_from=0)))
            region = TagNode(c.region, c.spec.id)
            if c.extent is not None:
                region = SpanIntervalNode(region, *c.extent)
        else:
            region = TagNode(c.region, c.spec.id)
            if isinstance(c, TransverseContribution):
                region = SpanIntervalNode(region, c.z - 0.5 * c.thickness, c.z + 0.5 * c.thickness)
        material.append(region)
        priority.append(c.spec.priority)

    return Structure(
        table=ctx.table,
        span=ctx.span,
        oml=oml,
        thickness=thickness,
        wall_mid=wall_mid,
        cavity=cavity,
        wall_meta=wall_meta,
        material=TrimNode(PriorityUnionNode(material, priority), oml),
        walls=walls,
        branches=tuple(branches),
        regions=tuple(c for c in contributions if isinstance(c, RegionContribution)),
        transverse=tuple(sorted((c for c in contributions if isinstance(c, TransverseContribution)),
                                key=lambda t: t.z)),
    )


__all__ = ["BranchField", "Structure", "StructureAssembler", "StructureError", "WallThickness", "assemble"]
