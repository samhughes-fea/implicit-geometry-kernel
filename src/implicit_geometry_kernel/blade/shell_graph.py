"""Loftable shell graph, mid-surface partitions and topology events.

Built from a ``SectionCut``: caps are partitions on the skin loop, webs are
branches.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Dict, Tuple

import numpy as np

if TYPE_CHECKING:
    from implicit_geometry_kernel.blade.cut import SectionCut, SectionWall


@dataclass(frozen=True)
class MidSegment:
    """One straight piece of a wall mid-surface: point indices, name and wall layers (inside first)."""

    start: int
    end: int
    name: str
    kind: str  # "skin" or "web" ("skin" also for a cap merged into its skin)
    layers: Tuple[str, ...]  # SectionWall labels, inside first: (cap, skin) along a bond


@dataclass(frozen=True)
class MidSurface:
    """The walls of a ``SectionCut`` as a graph of oriented mid-surface segments."""

    points: np.ndarray  # (n, 2) blade-frame x, y
    segments: Tuple[MidSegment, ...]
    z: float = 0.0

    @property
    def length(self) -> np.ndarray:
        P = self.points
        return np.array([float(np.linalg.norm(P[s.end] - P[s.start])) for s in self.segments])


@dataclass(frozen=True)
class ShellRun:
    """One loftable run of the shell graph (a partition on the wall, or a web).

    ``surface`` is ``skin``, ``web``, ``rib``, or opt-in ``core`` / ``insert``.
    ``kind`` is the colour / feature class: ``cap`` on a stacked skin partition,
    the member name on a rib face, otherwise ``surface``.
    """

    name: str
    surface: str
    kind: str
    layers: Tuple[str, ...]
    points: np.ndarray
    thickness: float
    requires_eccentricity: bool = False
    eccentricity_m: float = 0.0
    layer_thicknesses_m: Tuple[float, ...] = ()
    skin_side: str | None = None
    reference_surface: str = "skin_mid"


@dataclass(frozen=True)
class PartitionEdge:
    """Feature line between two adjacent loop partitions at one station.

    ``point`` is the meeting vertex in the blade frame ``(x, y, z)``. Cap edges
    are the pairs whose names include a stacked cap run.
    """

    z: float
    name_a: str
    name_b: str
    point: np.ndarray


@dataclass(frozen=True)
class ShellGraph:
    """Loftable shell: one skin loop of partitions, web branches, and rib faces.

    Rib faces sit in ``faces`` at the member's ``z`` only. They are not lofted
    with the skin loop.
    """

    z: float
    loop: Tuple[ShellRun, ...]
    webs: Tuple[ShellRun, ...]
    faces: Tuple[ShellRun, ...] = ()
    branches: Tuple[ShellRun, ...] = ()
    solids: Tuple[ShellRun, ...] = ()

    def runs(self) -> Tuple[ShellRun, ...]:
        return self.loop + self.webs + self.solids

    def partition_edges(self) -> Tuple[PartitionEdge, ...]:
        """``(z, name_a, name_b) → xy`` where two loop runs share an endpoint."""
        loop = self.loop
        if len(loop) < 2:
            return ()
        z = float(self.z)
        scale = max(1.0, max(float(np.max(np.abs(r.points[:, :2]))) for r in loop))
        tol = 1e-9 * scale
        out = []
        seen: set = set()
        for i, a in enumerate(loop):
            a_ends = (np.asarray(a.points[0], float), np.asarray(a.points[-1], float))
            for b in loop[i + 1:]:
                if a.name == b.name:
                    continue
                for pa in a_ends:
                    for pb in (np.asarray(b.points[0], float), np.asarray(b.points[-1], float)):
                        if float(np.linalg.norm(pa[:2] - pb[:2])) > tol:
                            continue
                        key = frozenset((a.name, b.name))
                        if key in seen:
                            continue
                        seen.add(key)
                        out.append(PartitionEdge(z, a.name, b.name, 0.5 * (pa + pb)))
        return tuple(out)

    def edge(self, name_a: str, name_b: str) -> PartitionEdge | None:
        want = {name_a, name_b}
        for e in self.partition_edges():
            if {e.name_a, e.name_b} == want:
                return e
        return None

    def closed_loop_chain(self, *, rel_tol: float = 1e-6):
        """Ordered ``(name, xy, reversed)`` ring, or ``None`` if partitions do not meet."""
        strips = [(r.name, np.asarray(r.points[:, :2], float)) for r in self.loop if len(r.points) >= 2]
        if not strips:
            return None
        scale = max(1.0, max(float(np.max(np.abs(xy))) for _, xy in strips))
        return _chain(strips, rel_tol * scale)

    def loop_polygon(self) -> np.ndarray:
        """Closed ``(N, 2)`` mid-loop from ordered partitions (the snap/query contour)."""
        strips = [(r.name, np.asarray(r.points[:, :2], float)) for r in self.loop if len(r.points) >= 2]
        if not strips:
            raise ValueError("shell graph has no loop runs")
        chain = self.closed_loop_chain(rel_tol=1e-9)
        if chain is None:
            ring = np.vstack([xy for _, xy in strips])
        else:
            ring = np.vstack([xy[::-1] if rev else xy for _, xy, rev in chain])
        if len(ring) and not np.allclose(ring[0], ring[-1]):
            ring = np.vstack([ring, ring[0]])
        return ring

    def with_web_feet_on_loop(self) -> "ShellGraph":
        """Shared nodes: web-midline ends lie on the skin-mid loop (not T-joint feet).

        2D-F already shares those vertices (gap 0). This is the flush contract
        on the graph, not a solids recipe change.
        """
        from dataclasses import replace

        if not self.loop or not self.webs:
            return self
        loop_xy = np.vstack([np.asarray(r.points, float) for r in self.loop])
        webs = []
        for web in self.webs:
            pts = np.array(web.points, float, copy=True)
            for i in (0, -1):
                j = int(np.argmin(np.linalg.norm(loop_xy[:, :2] - pts[i, :2], axis=1)))
                pts[i] = loop_xy[j]
            webs.append(replace(web, points=pts))
        return replace(self, webs=tuple(webs))

    def offset_regions(self) -> Tuple:
        """Loop partitions that need constitutive eccentricity (cap footprints)."""
        from implicit_geometry_kernel.contract.laminate_offset_region import LaminateOffsetRegion

        return tuple(
            LaminateOffsetRegion(
                name=run.name,
                requires_eccentricity=True,
                eccentricity_m=float(run.eccentricity_m),
                layers=run.layers,
                layer_thicknesses_m=run.layer_thicknesses_m,
                skin_side=run.skin_side,
                reference_surface=run.reference_surface,
            )
            for run in self.loop
            if run.requires_eccentricity
        )


@dataclass(frozen=True)
class TopologyEvent:
    """A run present on only one of two stations: no label lerp across the gap."""

    kind: str
    name: str
    side: str  # "lost" (in the first graph only) or "gained"


@dataclass(frozen=True)
class TopologyDiff:
    """Shared loft names plus events when two station signatures differ."""

    shared: Tuple[str, ...]
    events: Tuple[TopologyEvent, ...]


def shared_loft_names(names_a, names_b) -> Tuple[str, ...]:
    """Names to loft between stations: intersection, order of ``names_a``. No lerp."""
    present = set(names_b)
    return tuple(n for n in names_a if n in present)


def _event_kind(run: ShellRun, side: str) -> str:
    if run.surface == "rib":
        return "lost_rib" if side == "lost" else "gained_rib"
    if run.surface == "stringer":
        return "lost_stringer" if side == "lost" else "gained_stringer"
    if run.surface == "web" or run.kind == "web":
        return "dropped_web" if side == "lost" else "added_web"
    if run.kind == "cap":
        return "cap_runout" if side == "lost" else "cap_onset"
    return "bay_merge" if side == "lost" else "bay_split"


def _events_for(a_runs, b_runs) -> list:
    by_a = {r.name: r for r in a_runs}
    by_b = {r.name: r for r in b_runs}
    events = []
    for name, run in by_a.items():
        if name not in by_b:
            events.append(TopologyEvent(_event_kind(run, "lost"), name, "lost"))
    for name, run in by_b.items():
        if name not in by_a:
            events.append(TopologyEvent(_event_kind(run, "gained"), name, "gained"))
    return events


def topology_diff(a: ShellGraph, b: ShellGraph) -> TopologyDiff:
    """Compare two shell graphs. Loft ``shared``; do not invent missing strips.

    Rib faces and partial-span stringers are topology events when present on
    only one station; they are not lofted (``shared`` is skin and web names).
    """
    shared = shared_loft_names([r.name for r in a.runs()], [r.name for r in b.runs()])
    events = _events_for(a.runs(), b.runs())
    events.extend(_events_for(a.faces, b.faces))
    events.extend(_events_for(a.branches, b.branches))
    return TopologyDiff(shared, tuple(events))


def _side_from_name(name: str) -> str | None:
    key = name.lower()
    if "suction" in key:
        return "suction_surface"
    if "pressure" in key:
        return "pressure_surface"
    return None


def shell_strips_from_graph(cut: SectionCut):
    """Shell hand-off: one ``ShellMidlineStrip`` per graph run (caps as partitions).

    Does not emit a second cap surface. ``cut.walls`` still lists solid caps.
    The strip type is unchanged (Wave 3); ``build_shell_midline_strips`` is
    the station skin+web emit, not a second cap surface.
    """
    from implicit_geometry_kernel.contract.shell_midline_strip import ShellMidlineStrip

    strips = []
    for run in cut.shell_graph().runs():
        side = _side_from_name(run.name)
        strips.append(
            ShellMidlineStrip(
                label=run.name,
                kind=run.kind,
                midline_twisted_pitch_axis_blade_frame=np.asarray(run.points[:, :2], float),
                thickness_m=float(run.thickness),
                closed=False,
                surface=side if run.kind == "cap" else None,
                skin_side=run.skin_side or (side if run.kind == "skin" else None),
                reference_face=run.reference_surface,
                layers=run.layers,
                layer_thicknesses_m=run.layer_thicknesses_m,
                requires_eccentricity=bool(run.requires_eccentricity),
                eccentricity_m=float(run.eccentricity_m),
            )
        )
    return tuple(strips)


def _chain(strips, tol):
    """Order open polylines into one closed loop: [(label, xy, reversed)], or None."""
    if not strips:
        return None
    left = list(strips)
    label, xy = left.pop(0)
    out = [(label, xy, False)]
    end = xy[-1]
    while left:
        for k, (lb, p) in enumerate(left):
            if np.linalg.norm(p[0] - end) <= tol:
                out.append((lb, p, False))
                end = p[-1]
                break
            if np.linalg.norm(p[-1] - end) <= tol:
                out.append((lb, p, True))
                end = p[0]
                break
        else:
            return None
        left.pop(k)
    if np.linalg.norm(end - xy[0]) > tol:
        return None
    return out


from implicit_geometry_kernel.blade.shell_graph_build import (  # noqa: E402
    _mid_surface,
    _rib_faces,
    _runs_of,
    _solid_runs,
    _stringer_branches,
)


def shell_graph(cut: SectionCut, include_solids: bool = False) -> ShellGraph:
    """Loftable shell graph: skin partitions, web branches, rib faces and stringers.

    ``include_solids`` (default off) adds core and insert runs from polygon
    boundaries so loft can see them without changing default 2D-F run names.
    """
    from dataclasses import replace

    ms = cut.mid_surface()
    skins = tuple(s for s in ms.segments if s.kind == "skin")
    webs = tuple(s for s in ms.segments if s.kind == "web")
    graph = ShellGraph(float(cut.z), _runs_of(ms, skins, cut, "skin"),
                       _runs_of(ms, webs, cut, "web")).with_web_feet_on_loop()
    solids = _solid_runs(cut) if include_solids else ()
    return replace(graph, faces=_rib_faces(cut, graph), branches=_stringer_branches(cut), solids=solids)


def partition_edge_curve(geo, name_a: str, name_b: str, zs) -> np.ndarray:
    """Spanwise polyline of a named partition edge (graph vertices, not Newton).

    Stations where the pair is missing or the cut is infeasible are omitted
    (no label lerp).
    """
    from implicit_geometry_kernel.contract import InfeasibleSectionGeometry

    out = []
    for z in np.asarray(zs, float).ravel():
        try:
            edge = geo.section_at(float(z)).shell_graph().edge(name_a, name_b)
        except InfeasibleSectionGeometry:
            continue
        if edge is not None:
            out.append(edge.point)
    return np.vstack(out) if out else np.zeros((0, 3), dtype=float)


__all__ = [
    "MidSegment",
    "MidSurface",
    "PartitionEdge",
    "ShellGraph",
    "ShellRun",
    "TopologyDiff",
    "TopologyEvent",
    "partition_edge_curve",
    "shared_loft_names",
    "shell_graph",
    "shell_strips_from_graph",
    "topology_diff",
]

