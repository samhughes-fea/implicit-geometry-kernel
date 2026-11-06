"""Loftable shell graph, mid-surface partitions and topology events.

Built from a ``SectionCut``: caps are partitions on the skin loop, webs are
branches. ``blade.cut`` re-exports these names until 0.6.0.
"""

from __future__ import annotations

from dataclasses import dataclass
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

    ``surface`` is ``skin`` or ``web`` (geometry). ``kind`` is the colour / feature
    class: ``cap`` on a stacked skin partition, otherwise ``surface``.
    """

    name: str
    surface: str
    kind: str
    layers: Tuple[str, ...]
    points: np.ndarray
    thickness: float


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
    """Loftable shell: one skin loop of partitions, plus web branches. No extra cap loop."""

    z: float
    loop: Tuple[ShellRun, ...]
    webs: Tuple[ShellRun, ...]

    def runs(self) -> Tuple[ShellRun, ...]:
        return self.loop + self.webs

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

    def with_web_feet_on_loop(self) -> "ShellGraph":
        """Snap web endpoints onto the nearest stacked-loop vertex (shell re-host).

        2D-F feet already share those vertices (gap 0); this is the explicit
        Family D / flush contract on the graph, not a solids recipe change.
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
        return ShellGraph(self.z, self.loop, tuple(webs))


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
    if run.surface == "web" or run.kind == "web":
        return "dropped_web" if side == "lost" else "added_web"
    if run.kind == "cap":
        return "cap_runout" if side == "lost" else "cap_onset"
    return "bay_merge" if side == "lost" else "bay_split"


def topology_diff(a: ShellGraph, b: ShellGraph) -> TopologyDiff:
    """Compare two shell graphs. Loft ``shared``; do not invent missing strips."""
    by_a = {r.name: r for r in a.runs()}
    by_b = {r.name: r for r in b.runs()}
    shared = shared_loft_names(by_a, by_b)
    events = []
    for name, run in by_a.items():
        if name not in by_b:
            events.append(TopologyEvent(_event_kind(run, "lost"), name, "lost"))
    for name, run in by_b.items():
        if name not in by_a:
            events.append(TopologyEvent(_event_kind(run, "gained"), name, "gained"))
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
    The strip type is unchanged (Wave 3); only the producer set differs from
    ``build_shell_midline_strips``.
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
                skin_side=side if run.kind == "skin" else None,
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


def _mid_surface(cut: SectionCut) -> MidSurface:
    tol = 1e-9 * max(1.0, max(float(np.abs(w.xy).max()) for w in cut.walls))
    pts = []
    lookup: Dict[Tuple[int, int], list] = {}

    def index(p):
        key = (int(round(p[0] / (100 * tol))), int(round(p[1] / (100 * tol))))
        for dk in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1)):
            for k in lookup.get((key[0] + dk[0], key[1] + dk[1]), ()):
                if np.linalg.norm(pts[k] - p) <= tol:
                    return k
        pts.append(np.array(p, float))
        lookup.setdefault(key, []).append(len(pts) - 1)
        return len(pts) - 1

    skins, webs, caps = cut.walls_of("skin"), cut.walls_of("web"), cut.walls_of("cap")
    # G9: cores / TE inserts are not SectionWalls; CoreWall.to_midline is not implemented.
    # skins: orient every skin so that z x s points out of the section
    flip = {w.label: False for w in skins}
    loop = _chain([(w.label, w.xy) for w in skins], tol)
    if loop is not None:
        ring = np.vstack([xy[::-1] if rev else xy for _, xy, rev in loop])
        area = 0.5 * np.sum(ring[:, 0] * np.roll(ring[:, 1], -1) - np.roll(ring[:, 0], -1) * ring[:, 1])
        for label, _, rev in loop:
            flip[label] = rev != (area > 0.0)  # counter-clockwise: z x s points in
    # caps: which skin vertices carry which cap, and where the stack's mid-surface lies
    capped: Dict[Tuple[str, int], str] = {}
    moved: Dict[Tuple[str, int], np.ndarray] = {}
    skin_of = {w.label: w for w in skins}
    for cap in caps:
        for host_label, (ic, ih) in cap.bonds.items():
            host = skin_of.get(host_label)
            if host is None:
                continue
            if len(ih) < 2 or np.sum(np.linalg.norm(np.diff(host.xy[ih], axis=0), axis=1)) < cap.thickness:
                # the cap's edge reaching just past a web onto the next skin: whether that
                # sliver holds a vertex depends on sampling, so it is left to the skin
                continue
            f = cap.thickness / (cap.thickness + host.thickness)
            for a, b in zip(ic, ih):
                capped[(host_label, int(b))] = cap.label
                moved[(host_label, int(b))] = host.xy[b] + f * (cap.xy[a] - host.xy[b])
            # likewise a sliver of bare skin between the cap's edge and the end of its skin
            arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(host.xy, axis=0), axis=1))])
            lo, hi = int(np.min(ih)), int(np.max(ih))
            for edge, rest, gap in ((lo, range(0, lo), arc[lo]), (hi, range(hi + 1, len(host.xy)), arc[-1] - arc[hi])):
                if len(rest) and gap < cap.thickness:
                    step = moved[(host_label, edge)] - host.xy[edge]
                    for v in rest:
                        capped[(host_label, v)] = cap.label
                        moved[(host_label, v)] = host.xy[v] + step
    segments = []
    shift: Dict[int, np.ndarray] = {}
    for w in skins:
        ids = [index(p) for p in w.xy]
        for v, k in enumerate(ids):
            if (w.label, v) in moved:
                shift[k] = moved[(w.label, v)]
        run = []
        for v in range(len(ids) - 1):
            i, j = ids[v], ids[v + 1]
            if i == j:
                continue
            ca, cb = capped.get((w.label, v)), capped.get((w.label, v + 1))
            if ca is not None and ca == cb:
                seg = (i, j, ca, (ca, w.label))
            else:
                seg = (i, j, w.label, (w.label,))
            run.append(seg)
        if flip[w.label]:
            run = [(j, i, name, layers) for i, j, name, layers in reversed(run)]
        segments.extend(MidSegment(i, j, name, "skin", layers) for i, j, name, layers in run)
    t = np.radians(cut.twist)
    chord_dir = np.array([np.cos(t), -np.sin(t)])
    for w in webs:
        ids = [index(p) for p in w.xy]
        run = [(i, j) for i, j in zip(ids[:-1], ids[1:]) if i != j]
        if run:
            d = pts[run[-1][1]] - pts[run[0][0]]
            if np.array([-d[1], d[0]]) @ chord_dir < 0.0:
                run = [(j, i) for i, j in reversed(run)]
        segments.extend(MidSegment(i, j, w.label, "web", (w.label,)) for i, j in run)
    P = np.array(pts)
    for k, p in shift.items():
        P[k] = p
    return MidSurface(P, tuple(segments), cut.z)


def _colour_kind(surface: str, layers: Tuple[str, ...], walls: Dict[str, SectionWall]) -> str:
    if surface == "web":
        return "web"
    for label in layers:
        w = walls.get(label)
        if w is not None and w.kind == "cap":
            return "cap"
    return "skin"


def _runs_of(ms: MidSurface, segs, cut: SectionCut, surface: str) -> Tuple[ShellRun, ...]:
    if not segs:
        return ()
    by_label = {w.label: w for w in cut.walls}
    z = float(cut.z)

    def xyz(i: int) -> np.ndarray:
        p = ms.points[int(i)]
        return np.array([p[0], p[1], z], dtype=float)

    def pack(name, layers, pts) -> ShellRun:
        t = sum(by_label[lb].thickness for lb in layers if lb in by_label)
        return ShellRun(name, surface, _colour_kind(surface, layers, by_label), layers,
                        np.vstack(pts), float(t))

    runs = []
    name, layers = segs[0].name, segs[0].layers
    pts = [xyz(segs[0].start), xyz(segs[0].end)]
    for s in segs[1:]:
        if s.name == name and s.layers == layers:
            pts.append(xyz(s.end))
            continue
        runs.append(pack(name, layers, pts))
        name, layers = s.name, s.layers
        pts = [xyz(s.start), xyz(s.end)]
    runs.append(pack(name, layers, pts))
    return tuple(runs)


def shell_graph(cut: SectionCut) -> ShellGraph:
    """Loftable shell graph: skin partitions (caps merged) and web branches."""
    ms = cut.mid_surface()
    skins = tuple(s for s in ms.segments if s.kind == "skin")
    webs = tuple(s for s in ms.segments if s.kind == "web")
    return ShellGraph(float(cut.z), _runs_of(ms, skins, cut, "skin"),
                      _runs_of(ms, webs, cut, "web")).with_web_feet_on_loop()


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

