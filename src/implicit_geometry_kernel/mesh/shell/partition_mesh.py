"""The (t, z) partition mesher: a structured quad grid on every partition (decision 0015).

:func:`partition_mesh` meshes each skin and web partition of the shell graph as
its own (t, z) patch: ``t`` across the partition by normalised arc length
(:func:`~implicit_geometry_kernel.blade.skeleton.run_t`), ``z`` along the span
at the mesh stations. Each patch is a structured grid of quads, so element
edges follow the partition's own coordinates.

Nodes are shared by skeleton identity, not welded by coordinate: a partition
boundary is one node per station, used by both partitions, and a web foot is
the skin node at that boundary or vertex. Mesh lines run through every
skeleton point on a partition (web feet) and through every sub-partition band
edge (decision 0016), so a thickness or eccentricity step lies on an element
edge.

Nodes lie on the section mid-line at each station (the cut's run polylines).
This first step meshes blades whose shell topology is the same at every
station, with a closed skin loop, and builds skin and webs only. Topology
changes (caps or webs that start or stop inside the span, as transition
bands), the rib and tip-closure fill, and Newton junction nodes between
stations are later steps of plan row L6. It sits beside :func:`shell_mesh`,
which stays the default.
"""

from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np

from implicit_geometry_kernel.blade.skeleton import run_point_at
from implicit_geometry_kernel.mesh.shell.orient import _oriented
from implicit_geometry_kernel.mesh.shell.shell_mesh import ShellMesh
from implicit_geometry_kernel.mesh.shell.stations import mesh_stations

_T_TOL = 1e-9


class PartitionMeshError(ValueError):
    """The blade is outside what the partition mesher handles yet."""


def _graph_signature(graph) -> tuple:
    chain = graph.closed_loop_chain()
    if chain is None:
        raise PartitionMeshError(
            f"station z={float(graph.z):g}: the skin loop is not closed; the partition mesher needs a closed loop"
        )
    if any(rev for _, _, rev in chain):
        raise PartitionMeshError(
            f"station z={float(graph.z):g}: a loop partition runs against the loop direction"
        )
    if graph.branches:
        raise PartitionMeshError(
            f"station z={float(graph.z):g}: branches ({', '.join(b.name for b in graph.branches)}) are not meshed yet"
        )
    return tuple(name for name, _, _ in chain), tuple(sorted(w.name for w in graph.webs))


def _breakpoints(skeleton, run, bands: Sequence[Any]) -> dict:
    """``{key: t}`` of interior points on ``run`` that a mesh line must pass through."""
    out: dict = {}
    for node in skeleton.nodes:
        t = node.runs.get(run.name)
        if t is not None and _T_TOL < t < 1.0 - _T_TOL:
            out[node.key] = float(t)
    for s in bands:
        if s.host != run.name:
            continue
        for k, t in enumerate(s.t_range):
            if _T_TOL < t < 1.0 - _T_TOL:
                out[("band", s.name, k)] = float(t)
    return out


def _segment_breaks(t_by_key: dict) -> list:
    """Sorted ``[(key, t)]`` with the run ends; coincident breakpoints merged."""
    items = sorted(t_by_key.items(), key=lambda kv: kv[1])
    merged: list = [(("end", 0), 0.0)]
    for key, t in items:
        if abs(t - merged[-1][1]) <= _T_TOL:
            continue
        merged.append((key, t))
    if abs(1.0 - merged[-1][1]) <= _T_TOL:
        merged.pop()
    merged.append((("end", 1), 1.0))
    return merged


def _te_point(ring: np.ndarray) -> np.ndarray:
    """The loop point farthest from the loop's mean point (the trailing edge of an aerofoil)."""
    pts = ring[:-1] if np.allclose(ring[0], ring[-1]) else ring
    c = pts.mean(axis=0)
    return pts[int(np.argmax(np.linalg.norm(pts - c, axis=1)))]


def partition_mesh(geo, *, h: float, zs: Sequence[float] | None = None,
                   sub_partitions: Sequence[Any] = ()) -> ShellMesh:
    """A structured (t, z) quad mesh of the blade's skin and webs.

    Parameters
    ----------
    geo:
        The blade (anything with ``section_at(z).shell_graph()`` and stations).
    h:
        Target element edge in metres, across each partition. Each segment
        between mesh lines gets the largest count it needs at any station, so
        the edge is at most ``h`` everywhere.
    zs:
        Span stations; default :func:`mesh_stations` with ``h``.
    sub_partitions:
        :class:`~implicit_geometry_kernel.structure.SubPartition` bands whose
        ``t`` edges get mesh lines. Their span windows are not mesh lines.

    Raises
    ------
    PartitionMeshError
        If the shell topology changes along the span or the skin loop is open.
    """
    if not float(h) > 0.0:
        raise ValueError("element size h must be positive")
    z_list = [float(z) for z in (zs if zs is not None else mesh_stations(geo, h=float(h)))]
    if len(z_list) < 2:
        raise ValueError("the partition mesher needs at least two stations")
    graphs = [geo.section_at(z).shell_graph() for z in z_list]
    signature = _graph_signature(graphs[0])
    for g in graphs[1:]:
        if _graph_signature(g) != signature:
            raise PartitionMeshError(
                f"shell topology changes between z={z_list[0]:g} and z={float(g.z):g}; "
                "the partition mesher handles constant topology for now (use shell_mesh)"
            )
    loop_names, web_names = signature
    hosts = {s.host for s in sub_partitions}
    unknown = hosts - set(loop_names) - set(web_names)
    if unknown:
        raise ValueError(f"sub-partition hosts not in the shell graph: {sorted(unknown)}")

    def runs_of(g) -> dict:
        return {r.name: r for r in g.loop + g.webs}

    by_station = [runs_of(g) for g in graphs]
    skeletons = [g.skeleton() for g in graphs]
    order = list(loop_names) + list(web_names)

    # Mesh lines on each run: the same keys, in the same order, at every station.
    breaks: dict = {}
    counts: dict = {}
    for name in order:
        per_station = [_segment_breaks(_breakpoints(sk, rs[name], sub_partitions)) for sk, rs in zip(skeletons, by_station)]
        keys0 = [k for k, _ in per_station[0]]
        for g, b in zip(graphs[1:], per_station[1:]):
            if [k for k, _ in b] != keys0:
                raise PartitionMeshError(
                    f"partition {name!r}: its mesh lines change order or number at z={float(g.z):g}"
                )
        n_seg = len(keys0) - 1
        need = [1] * n_seg
        for rs, b in zip(by_station, per_station):
            pts = np.asarray(rs[name].points, float)
            length = float(np.sum(np.linalg.norm(np.diff(pts[:, :2], axis=0), axis=1)))
            for k in range(n_seg):
                need[k] = max(need[k], int(math.ceil((b[k + 1][1] - b[k][1]) * length / float(h) - 1e-9)))
        breaks[name] = per_station
        counts[name] = need

    nodes: list[np.ndarray] = []

    def new_node(p: np.ndarray) -> int:
        nodes.append(np.asarray(p, float))
        return len(nodes) - 1

    ids: dict = {}  # (station, run) -> node ids across the run
    for j, (g, rs) in enumerate(zip(graphs, by_station)):
        z = z_list[j]
        line_ids: dict = {}  # skeleton key -> node id at this station
        prev_end: int | None = None
        first_start: int | None = None
        for name in loop_names:
            run = rs[name]
            pts = np.asarray(run.points, float)
            b = breaks[name][j]
            row: list[int] = []
            closing = name == loop_names[-1] and first_start is not None
            for k in range(len(b) - 1):
                ts = np.linspace(b[k][1], b[k + 1][1], counts[name][k] + 1)
                for m, t in enumerate(ts[:-1] if k < len(b) - 2 else ts):
                    if k == 0 and m == 0 and prev_end is not None:
                        row.append(prev_end)
                        continue
                    if closing and first_start is not None and k == len(b) - 2 and m == len(ts) - 1:
                        row.append(first_start)  # the ring closes on the first run's start
                        continue
                    key = b[k][0] if m == 0 else None
                    p = run_point_at(pts, float(t))
                    nid = new_node(np.array([p[0], p[1], z]))
                    row.append(nid)
                    if key is not None:
                        line_ids[key] = nid
            if first_start is None:
                first_start = row[0]
            prev_end = row[-1]
            ids[(j, name)] = row
        for name in web_names:
            run = rs[name]
            pts = np.asarray(run.points, float)
            b = breaks[name][j]
            row = []
            ends = {}
            sk = skeletons[j].by_key()
            for end in (0, 1):
                node = sk.get(("junction", name, end))
                if node is None:
                    continue
                for host in loop_names:
                    t = node.runs.get(host)
                    if t is None:
                        continue
                    hrow = ids[(j, host)]
                    if t <= _T_TOL:
                        ends[end] = hrow[0]
                    elif t >= 1.0 - _T_TOL:
                        ends[end] = hrow[-1]
                    elif node.key in line_ids:
                        ends[end] = line_ids[node.key]
                    if end in ends:
                        break
            for k in range(len(b) - 1):
                ts = np.linspace(b[k][1], b[k + 1][1], counts[name][k] + 1)
                for m, t in enumerate(ts[:-1] if k < len(b) - 2 else ts):
                    first = k == 0 and m == 0
                    last_pt = k == len(b) - 2 and m == len(ts) - 1
                    if first and 0 in ends:
                        row.append(ends[0])
                    elif last_pt and 1 in ends:
                        row.append(ends[1])
                    else:
                        p = run_point_at(pts, float(t))
                        row.append(new_node(np.array([p[0], p[1], z])))
            ids[(j, name)] = row

    elements: list[tuple[int, ...]] = []
    rows: list[tuple] = []
    for j in range(len(z_list) - 1):
        for name in order:
            a, b = ids[(j, name)], ids[(j + 1, name)]
            r0, r1 = by_station[j][name], by_station[j + 1][name]
            th = 0.5 * (float(r0.thickness) + float(r1.thickness))
            ecc = 0.5 * (float(r0.eccentricity_m) + float(r1.eccentricity_m))
            surface = "web" if name in web_names else "skin"
            for i in range(len(a) - 1):
                elements.append((a[i], a[i + 1], b[i + 1], b[i]))
                rows.append((name, surface, th, ecc, tuple(r0.layers), tuple(r0.layer_thicknesses_m), j,
                             z_list[j], z_list[j + 1]))

    pts = np.asarray(nodes, float)
    te_at = {z: np.append(_te_point(np.asarray(g.loop_polygon(), float)), z) for z, g in zip(z_list, graphs)}
    elements, flipped = _oriented(pts, elements, [r[1] for r in rows], z_list, te_at)
    return ShellMesh(
        nodes=pts,
        elements=tuple(tuple(int(i) for i in e) for e in elements),
        name=tuple(r[0] for r in rows),
        surface=tuple(r[1] for r in rows),
        partition=tuple(r[0] for r in rows),
        thickness=np.array([r[2] for r in rows], float),
        eccentricity=np.array([r[3] for r in rows], float),
        layers=tuple(r[4] for r in rows),
        layer_thicknesses=tuple(r[5] for r in rows),
        stations=tuple(z_list),
        n_source=len(pts),
        bay=tuple(r[6] for r in rows),
        z0=tuple(r[7] for r in rows),
        z1=tuple(r[8] for r in rows),
        constructor=tuple("partition" for _ in rows),
        flipped=tuple(int(f) for f in flipped),
        node_origin=tuple("station" for _ in range(len(pts))),
        section_distance=np.zeros(len(pts)),
    )


__all__ = ["PartitionMeshError", "partition_mesh"]
