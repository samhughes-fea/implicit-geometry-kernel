"""The shell skeleton: tagged curves where partitions and members meet (decision 0015).

At one station, :func:`section_skeleton` turns a :class:`~implicit_geometry_kernel.blade.shell_graph.ShellGraph`
into :class:`SkeletonNode` objects, one per point where runs meet or end. Each node
has one or more classes:

``junction``
    A web end on the skin loop.
``partition_boundary``
    Two loop partitions meet (cap edges, skin splits at the leading and
    trailing edge).
``free_edge``
    A run end that meets nothing (a stringer's free end, an open trailing edge).

Along the span, :func:`blade_skeleton` follows each node by its key through a
list of stations and returns :class:`SkeletonCurve` objects. A curve that is not
present at the first or last station starts or ends at a ``termination`` (a cap
or web that starts or stops inside the span).

Each run also gets a normalised coordinate ``t`` in [0, 1] by arc length along
its points (:func:`run_t`), the across coordinate of the run's (t, z) patch.
Junction nodes carry the in-plane web-to-skin angle ``phi_deg`` (90 for a web
normal to the skin), measured from the skin tangent in the loop's own direction
to the web direction pointing away from the foot.

Geometry only: names, points in metres, angles in degrees.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence, Tuple

import numpy as np

SKELETON_CLASSES = ("junction", "partition_boundary", "free_edge", "termination")

#: Web-to-skin angles outside this range (degrees) are flagged by ``validate()``.
PHI_RANGE_DEG = (30.0, 150.0)


def run_t(points: np.ndarray) -> np.ndarray:
    """Normalised arc length in [0, 1] at each point of a run (0 at the first point).

    A run of zero length gives zeros.
    """
    pts = np.asarray(points, float)
    if len(pts) < 2:
        return np.zeros(len(pts))
    seg = np.linalg.norm(np.diff(pts[:, :2], axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(s[-1])
    return s / total if total > 0.0 else np.zeros(len(pts))


def run_point_at(points: np.ndarray, t: float) -> np.ndarray:
    """Point at normalised arc length ``t`` along a run (linear between its points)."""
    pts = np.asarray(points, float)
    tt = run_t(pts)
    tc = float(np.clip(t, 0.0, 1.0))
    return np.array([np.interp(tc, tt, pts[:, k]) for k in range(pts.shape[1])])


@dataclass(frozen=True)
class SkeletonNode:
    """One point of the skeleton at one station.

    Attributes
    ----------
    key:
        Stable identity along the span: ``("junction", web, end)`` with ``end``
        0 or 1, ``("partition_boundary", a, b)`` with ``a < b``, or
        ``("free_edge", run, end)``.
    classes:
        The node's classes (a web foot at a cap edge is both a junction and a
        partition boundary).
    point:
        ``(x, y, z)`` in the blade frame, metres.
    runs:
        ``{run name: t}`` for every run that has this point as an end or vertex.
    phi_deg:
        In-plane web-to-skin angle in degrees at a junction, else ``None``.
    """

    key: Tuple
    classes: frozenset
    point: np.ndarray
    runs: Mapping[str, float] = field(default_factory=dict)
    phi_deg: float | None = None


@dataclass(frozen=True)
class SectionSkeleton:
    """The skeleton nodes of one station."""

    z: float
    nodes: Tuple[SkeletonNode, ...]

    def of_class(self, name: str) -> Tuple[SkeletonNode, ...]:
        """Nodes carrying class ``name``."""
        return tuple(n for n in self.nodes if name in n.classes)

    def by_key(self) -> dict:
        """``{key: node}``."""
        return {n.key: n for n in self.nodes}


@dataclass(frozen=True)
class SkeletonCurve:
    """One skeleton node followed along the span.

    ``z`` and ``points`` hold the stations where the node exists. ``starts`` and
    ``ends`` say whether each end is a ``termination`` (the node appears or
    disappears inside the span) or the blade's first or last station.
    """

    key: Tuple
    classes: frozenset
    z: np.ndarray
    points: np.ndarray
    phi_deg: np.ndarray
    starts: str
    ends: str


def _scale(graph) -> float:
    pts = [np.asarray(r.points, float) for r in graph.loop + graph.webs + graph.branches if len(r.points)]
    if not pts:
        return 1.0
    return max(1.0, float(max(np.max(np.abs(p[:, :2])) for p in pts)))


def _loop_ring(graph) -> np.ndarray | None:
    if not graph.loop:
        return None
    try:
        ring = graph.loop_polygon()
    except ValueError:
        return None
    return np.asarray(ring, float)


def _tangent_at(ring: np.ndarray, xy: np.ndarray, tol: float) -> np.ndarray | None:
    pts = ring[:-1] if np.allclose(ring[0], ring[-1]) else ring
    d = np.linalg.norm(pts[:, :2] - xy[:2], axis=1)
    i = int(np.argmin(d))
    if float(d[i]) > tol:
        return None
    n = len(pts)
    t = pts[(i + 1) % n, :2] - pts[(i - 1) % n, :2]
    norm = float(np.linalg.norm(t))
    return t / norm if norm > 0.0 else None


def _phi_deg(ring, web_pts: np.ndarray, end: int, tol: float) -> float | None:
    if ring is None or len(web_pts) < 2:
        return None
    foot = web_pts[0] if end == 0 else web_pts[-1]
    inward = (web_pts[1] if end == 0 else web_pts[-2])[:2] - foot[:2]
    norm = float(np.linalg.norm(inward))
    tangent = _tangent_at(ring, foot, tol)
    if tangent is None or norm == 0.0:
        return None
    cos = float(np.clip(np.dot(tangent, inward / norm), -1.0, 1.0))
    return float(np.degrees(np.arccos(cos)))


def _runs_at(graph, xy: np.ndarray, tol: float) -> dict:
    out: dict = {}
    for run in graph.loop + graph.webs + graph.branches:
        pts = np.asarray(run.points, float)
        if len(pts) == 0:
            continue
        d = np.linalg.norm(pts[:, :2] - xy[:2], axis=1)
        i = int(np.argmin(d))
        if float(d[i]) <= tol:
            out[run.name] = float(run_t(pts)[i])
    return out


def section_skeleton(graph) -> SectionSkeleton:
    """Tagged skeleton nodes of one :class:`ShellGraph` (one station)."""
    scale = _scale(graph)
    tol = 1e-9 * scale
    z = float(graph.z)
    ring = _loop_ring(graph)
    nodes: dict = {}

    def add(key, cls, point):
        p = np.asarray(point, float)
        if len(p) == 2:
            p = np.array([p[0], p[1], z])
        node = nodes.get(key)
        if node is None:
            nodes[key] = SkeletonNode(key, frozenset({cls}), p, _runs_at(graph, p, tol))
        else:
            nodes[key] = SkeletonNode(key, node.classes | {cls}, node.point, node.runs, node.phi_deg)

    for edge in graph.partition_edges():
        a, b = sorted((edge.name_a, edge.name_b))
        add(("partition_boundary", a, b), "partition_boundary", edge.point)

    for web in graph.webs:
        pts = np.asarray(web.points, float)
        if len(pts) < 2:
            continue
        for end in (0, 1):
            foot = pts[0] if end == 0 else pts[-1]
            key = ("junction", web.name, end)
            on_loop = ring is not None and _tangent_at(ring, foot, tol) is not None
            if not on_loop:
                add(("free_edge", web.name, end), "free_edge", foot)
                continue
            add(key, "junction", foot)
            node = nodes[key]
            nodes[key] = SkeletonNode(
                key, node.classes, node.point, node.runs, _phi_deg(ring, pts, end, tol)
            )

    for run in graph.branches:
        pts = np.asarray(run.points, float)
        if len(pts) < 2:
            continue
        for end in (0, 1):
            p = pts[0] if end == 0 else pts[-1]
            if ring is not None and _tangent_at(ring, p, tol) is not None:
                add(("junction", run.name, end), "junction", p)
            else:
                add(("free_edge", run.name, end), "free_edge", p)

    if graph.loop and graph.closed_loop_chain() is None:
        ends = {}
        for run in graph.loop:
            pts = np.asarray(run.points, float)
            for end in (0, 1):
                p = pts[0] if end == 0 else pts[-1]
                ends[(run.name, end)] = p
        for (name, end), p in ends.items():
            touching = sum(1 for (n2, _), q in ends.items() if n2 != name and np.linalg.norm(q[:2] - p[:2]) <= tol)
            if touching == 0:
                add(("free_edge", name, end), "free_edge", p)

    def merge_coincident(items):
        """A web foot at a cap edge is one point: carry both classes on each key."""
        out = list(items)
        for i, a in enumerate(out):
            for j, b in enumerate(out):
                if i != j and np.linalg.norm(a.point[:2] - b.point[:2]) <= tol:
                    phi = a.phi_deg if a.phi_deg is not None else b.phi_deg
                    out[i] = SkeletonNode(a.key, a.classes | b.classes, a.point, a.runs, phi)
                    a = out[i]
        return out

    return SectionSkeleton(z, tuple(merge_coincident(nodes.values())))


def blade_skeleton(geo, zs: Iterable[float]) -> Tuple[SkeletonCurve, ...]:
    """Skeleton curves of a blade, followed through stations ``zs`` (increasing).

    A node present at some stations and not at the first or last one starts or
    ends at a ``termination``; the curve then also carries that class.
    """
    z_list = [float(z) for z in zs]
    if not z_list:
        return ()
    skeletons = [section_skeleton(geo.section_at(z).shell_graph()) for z in z_list]
    keys: list = []
    for sk in skeletons:
        for node in sk.nodes:
            if node.key not in keys:
                keys.append(node.key)
    curves = []
    for key in keys:
        present = [(i, sk.by_key()[key]) for i, sk in enumerate(skeletons) if key in sk.by_key()]
        idx = [i for i, _ in present]
        classes = frozenset().union(*(n.classes for _, n in present))
        starts = "termination" if idx[0] > 0 else "blade_end"
        ends = "termination" if idx[-1] < len(z_list) - 1 else "blade_end"
        if "termination" in (starts, ends):
            classes = classes | {"termination"}
        curves.append(
            SkeletonCurve(
                key=key,
                classes=classes,
                z=np.array([z_list[i] for i in idx]),
                points=np.array([n.point for _, n in present]),
                phi_deg=np.array([np.nan if n.phi_deg is None else n.phi_deg for _, n in present]),
                starts=starts,
                ends=ends,
            )
        )
    return tuple(curves)


def phi_issues(skeleton: SectionSkeleton, phi_range: Sequence[float] = PHI_RANGE_DEG) -> list[str]:
    """Junctions whose web-to-skin angle lies outside ``phi_range`` (degrees)."""
    lo, hi = float(phi_range[0]), float(phi_range[1])
    out = []
    for node in skeleton.of_class("junction"):
        if node.phi_deg is not None and not lo <= node.phi_deg <= hi:
            out.append(
                f"station z={skeleton.z:g} junction {node.key[1]!r} end {node.key[2]}: "
                f"web-to-skin angle {node.phi_deg:.1f} deg outside [{lo:g}, {hi:g}]"
            )
    return out


__all__ = [
    "PHI_RANGE_DEG",
    "SKELETON_CLASSES",
    "SectionSkeleton",
    "SkeletonCurve",
    "SkeletonNode",
    "blade_skeleton",
    "phi_issues",
    "run_point_at",
    "run_t",
    "section_skeleton",
]
