"""Station wall network: labelled midlines with integer feature ids.

This is the geometric graph a mesher lofts. It is extracted from a ``SectionCut``,
not from Gmsh. Feature ids come from :class:`~implicit_geometry_kernel.sdf.nodes.base.FeatureTable`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from implicit_geometry_kernel.blade.cut import SectionCut, SectionWall
from implicit_geometry_kernel.sdf.nodes.base import NO_FEATURE, FeatureSpec, FeatureTable, SDFNode
from implicit_geometry_kernel.sdf.nodes.tagged import PolylineBandSDF, PriorityUnionNode, TagNode


@dataclass(frozen=True)
class WallBranch:
    """One wall of the station network (a loftable midline)."""

    spec: FeatureSpec
    kind: str
    points: np.ndarray
    thickness: float
    side: str | None = None

    @property
    def label(self) -> str:
        return self.spec.name

    @property
    def feature_id(self) -> int:
        return self.spec.id


@dataclass(frozen=True)
class StationNetwork:
    """Wall loop / branches at one ``z``, with ``loop_features`` for a mesher."""

    z: float
    walls: tuple[WallBranch, ...]

    @property
    def loop_features(self) -> tuple[int, ...]:
        return tuple(w.feature_id for w in self.walls)

    def wall_meta_field(self) -> SDFNode:
        """Priority-union of wall-strip SDFs tagged with each branch's feature id."""
        nodes: list[SDFNode] = []
        prios: list[int] = []
        for w in self.walls:
            band = PolylineBandSDF(w.points[:, :2], 0.5 * w.thickness, w.feature_id)
            nodes.append(TagNode(band, w.feature_id))
            prios.append(w.spec.priority)
        return PriorityUnionNode(nodes, prios)


def station_network(cut: SectionCut, table: FeatureTable | None = None) -> StationNetwork:
    """Build a station network from the merged shell graph (caps are partitions, not extra loops)."""
    from implicit_geometry_kernel.blade.shell_graph import shell_graph

    if table is None:
        table = FeatureTable()
    walls = []
    for run in shell_graph(cut).runs():
        spec = table.add(run.name, run.kind)
        walls.append(WallBranch(spec, run.kind, np.asarray(run.points, float),
                                float(run.thickness), None))
    return StationNetwork(float(cut.z), tuple(walls))


def wall_meta_from_cut(cut: SectionCut, table: FeatureTable | None = None) -> tuple[SDFNode, FeatureTable]:
    """Tagged bands of the shell graph: one wall loop plus webs, no nested cap strip."""
    from implicit_geometry_kernel.blade.shell_graph import shell_graph

    if table is None:
        table = FeatureTable()
    nodes: list[SDFNode] = []
    prios: list[int] = []
    for run in shell_graph(cut).runs():
        spec = table.add(run.name, run.kind)
        band = PolylineBandSDF(run.points[:, :2], max(0.5 * float(run.thickness), 1e-6), spec.id)
        nodes.append(TagNode(band, spec.id))
        prios.append(spec.priority)
    if not nodes:
        raise ValueError("section has no shell runs to tag")
    return PriorityUnionNode(nodes, prios), table


def station_network_from_meta(wall_mid: SDFNode, wall_meta: SDFNode, z: float,
                              bounds: tuple[float, float, float, float], table: FeatureTable,
                              n: int = 48) -> StationNetwork:
    """Station network v2: extract ``wall_mid`` loop, split by ``wall_meta`` ids.

    Thickness is left ``0.0``: this is a meta-split of the SDF loop, not the
    shell-graph gauge. Use ``station_network(cut)`` / ``ShellRun.thickness``.
    """
    from implicit_geometry_kernel.extract.marching import extract_zero_loops

    loops = extract_zero_loops(wall_mid, z, bounds, n=n)
    if not loops:
        raise ValueError(f"no wall-mid loop at z={z:g}")
    loop = loops[0]
    ids = np.asarray(wall_meta.evaluate(loop).meta, dtype=int)
    walls = []
    start = 0
    for i in range(1, len(ids) + 1):
        done = i == len(ids)
        if not done and ids[i] == ids[start]:
            continue
        fid = int(ids[start])
        if fid != NO_FEATURE and 0 <= fid < len(table):
            spec = table[fid]
            walls.append(WallBranch(spec, spec.kind, loop[start:i], 0.0, None))
        start = i
    if not walls:
        raise ValueError(f"wall_meta produced no tagged branches at z={z:g}")
    return StationNetwork(float(z), tuple(walls))


def feature_table_from_walls(walls: tuple[SectionWall, ...] | list[SectionWall]) -> FeatureTable:
    table = FeatureTable()
    for wall in walls:
        table.add(wall.label, wall.kind)
    return table
