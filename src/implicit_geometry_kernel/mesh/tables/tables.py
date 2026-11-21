"""Mesh tables on a shell mesh: edges, open loops, the failure query, and the CSV dump.

Quality stays in :data:`QUALITY_DTYPE`. These tables are the mesh itself: node ids,
coordinates, connectivity, and the decisions that built each face.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from implicit_geometry_kernel.mesh.quality import _corner_angles

ORIGINS = ("station", "span_mid", "centroid", "station_as_midpoint", "interior")
CONSTRUCTORS = ("ruled", "midspan_split", "step_patch", "transition", "fan", "bowtie",
                "omitted", "inherited", "interior")
KINDS = ("station", "spanwise", "interior")
MEETS = ("skin", "skin_web", "rib", "tip", "open", "other")

EDGE_DTYPE = np.dtype([
    ("a", np.int64), ("b", np.int64), ("use", np.int64),
    ("e0", np.int64), ("e1", np.int64), ("e2", np.int64),
    ("length", np.float64), ("dx", np.float64), ("dy", np.float64), ("dz", np.float64),
    ("dihedral_deg", np.float64), ("chain_yaw_deg", np.float64), ("chain_pitch_deg", np.float64),
    ("kind", "U16"), ("meet", "U16"), ("sagitta", np.float64),
])
"""One undirected edge. ``e0``..``e2`` are the element ids, ``-1`` when unused."""

CORNER_DTYPE = np.dtype([
    ("element_id", np.int64), ("corner", np.int64), ("node_id", np.int64),
    ("x", np.float64), ("y", np.float64), ("z", np.float64),
    ("angle_deg", np.float64), ("yaw_deg", np.float64), ("pitch_deg", np.float64), ("roll_deg", np.float64),
    ("dx_in", np.float64), ("dy_in", np.float64), ("dz_in", np.float64),
    ("dx_out", np.float64), ("dy_out", np.float64), ("dz_out", np.float64),
    ("prev_id", np.int64), ("next_id", np.int64),
    ("constructor", "U32"), ("bay", np.int64), ("region", "U64"), ("surface", "U16"),
])
"""One element corner. Columns match :class:`Corner`."""


def _edge_users(row) -> list[int]:
    return [int(row[name]) for name in ("e0", "e1", "e2") if int(row[name]) >= 0]


@dataclass(frozen=True)
class StationDecision:
    """Why a span station is in the mesh, or why a grade refused to add one."""

    z: float
    origin: str
    refused: bool
    detail: str


@dataclass(frozen=True)
class StripDecision:
    """What one lofted strip did in one bay, including which strip made the bay an event."""

    bay: int
    z0: float
    z1: float
    region: str
    n0: int
    n1: int
    constructor: str
    trigger: str
    adv_a: int
    adv_b: int
    ties: int
    note: str
    chord_ratio: float = 0.0
    """``|c1 - c0| / max(c0, c1)`` of the bay. Zero when the chord does not change."""
    site0: int = -1
    """Coarse-edge index of the first 2:1 site. ``-1`` when the strip has no site."""
    site1: int = -1
    """Coarse-edge index of the second 2:1 site. ``-1`` when the strip has no site."""


@dataclass(frozen=True)
class MeshEdge:
    """One undirected pair of node ids."""

    a: int
    b: int
    use: int
    elements: tuple[int, ...]
    length: float
    dx: float
    dy: float
    dz: float
    """Edge vector from node ``a`` to node ``b``, in metres. ``a`` is the lower node id."""
    dihedral_deg: float
    """Signed fold from one face normal to the other about the edge vector, degrees.

    The pair is the two faces on the same surface when the edge meets more than
    two, otherwise the first two faces. The sign runs from the lower element id
    to the higher, about the vector from node ``a`` to node ``b``. Zero on a
    boundary edge and on a coplanar pair. The included angle where two skins
    meet is a steady value of this fold. A crinkle is ``chain_pitch_deg``.
    """
    chain_yaw_deg: float
    """Turn in the face plane from the previous spanwise edge of the same regions, degrees.

    Zero when the edge is not spanwise or is the first edge of its chain.
    """
    chain_pitch_deg: float
    """Turn out of the face plane from the previous spanwise edge, degrees.

    A trailing-edge crinkle is a jump in this angle from one spanwise edge to
    the next. Zero on the first edge of a chain.
    """
    kind: str
    meet: str
    sagitta: float


@dataclass(frozen=True)
class OpenLoop:
    """A chain of edges used by one element. The root ring is the expected one."""

    z: float
    n_edges: int
    is_root: bool
    nodes: tuple[int, ...]


@dataclass(frozen=True)
class Corner:
    """One interior angle of one element, named by the node it sits on.

    A node shared by four quads has four angles, so the angle is not a node
    attribute. ``x``, ``y``, ``z`` are that node's coordinates, copied onto the
    row so the corner can be plotted without a join.
    """

    element_id: int
    corner: int
    node_id: int
    x: float
    y: float
    z: float
    angle_deg: float
    """Interior angle in the plane of the two edges, degrees. This is the in-face corner."""
    yaw_deg: float
    """Turn of the outgoing edge in the element plane, degrees. Positive is a left turn about the normal."""
    pitch_deg: float
    """Turn of the outgoing edge out of the element plane, degrees. A crinkle lifts this off zero."""
    roll_deg: float
    """Twist of the surface normal about the outgoing edge, degrees, onto the face that shares that edge."""
    dx_in: float
    dy_in: float
    dz_in: float
    """Incoming edge vector, previous node to this node, metres."""
    dx_out: float
    dy_out: float
    dz_out: float
    """Outgoing edge vector, this node to the next, metres."""
    prev_id: int
    next_id: int
    constructor: str
    bay: int
    region: str
    surface: str


@dataclass(frozen=True)
class Seam:
    """Two node ids that did not merge, or a node lying in the interior of an edge."""

    kind: str
    a: int
    b: int
    node: int
    distance: float


def _meet(surfaces: set[str], use: int) -> str:
    if use <= 1:
        return "open"
    if "tip_closure" in surfaces:
        return "tip"
    if "rib" in surfaces:
        return "rib"
    if "skin" in surfaces and "web" in surfaces:
        return "skin_web"
    if surfaces == {"skin"}:
        return "skin"
    return "other"


def _kind(origin_a: str, origin_b: str, za: float, zb: float) -> str:
    if origin_a == "station" and origin_b == "station":
        return "station" if abs(za - zb) <= 1e-9 else "spanwise"
    return "interior"


def _unit(vector: np.ndarray):
    length = float(np.linalg.norm(vector))
    if length <= 1e-15:
        return None
    return vector / length


def _face_normal(points: np.ndarray) -> np.ndarray:
    """Unit Newell normal of an element loop. Zero when the loop has no area."""
    normal = np.zeros(3)
    for i, point in enumerate(points):
        nxt = points[(i + 1) % len(points)]
        normal[0] += (point[1] - nxt[1]) * (point[2] + nxt[2])
        normal[1] += (point[2] - nxt[2]) * (point[0] + nxt[0])
        normal[2] += (point[0] - nxt[0]) * (point[1] + nxt[1])
    unit = _unit(normal)
    return normal if unit is None else unit


def _trajectory_angles(incoming: np.ndarray, outgoing: np.ndarray,
                       normal_here: np.ndarray, normal_next: np.ndarray) -> tuple[float, float, float]:
    """``(yaw, pitch, roll)`` in degrees that carry the edge frame onto the next edge.

    Yaw turns the edge in the element plane, pitch turns it out of that plane,
    and roll twists the normal about the outgoing edge onto the next face.
    """
    yaw, pitch, roll = _trajectory_block(
        np.asarray(incoming, float).reshape(1, 3), np.asarray(outgoing, float).reshape(1, 3),
        np.asarray(normal_here, float).reshape(1, 3), np.asarray(normal_next, float).reshape(1, 3),
    )
    return float(yaw[0]), float(pitch[0]), float(roll[0])


def _unit_rows(vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    length = np.linalg.norm(vectors, axis=1)
    ok = length > 1e-15
    unit = np.divide(vectors, length[:, None], out=np.zeros_like(vectors), where=ok[:, None])
    return unit, ok


def _trajectory_block(incoming: np.ndarray, outgoing: np.ndarray,
                      normal_here: np.ndarray, normal_next: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``_trajectory_angles`` for one row per sample. A missing frame stays at zero."""
    tangent_here, ok_here = _unit_rows(incoming)
    tangent_there, ok_there = _unit_rows(outgoing)
    lifted_here = normal_here - tangent_here * np.sum(normal_here * tangent_here, axis=1)[:, None]
    lifted_there = normal_next - tangent_there * np.sum(normal_next * tangent_there, axis=1)[:, None]
    normal_u_here, ok_normal_here = _unit_rows(lifted_here)
    normal_u_there, ok_normal_there = _unit_rows(lifted_there)
    binormal_here = np.cross(normal_u_here, tangent_here)
    binormal_there = np.cross(normal_u_there, tangent_there)
    along = np.sum(normal_u_here * tangent_there, axis=1)
    yaw = np.degrees(np.arctan2(np.sum(binormal_here * tangent_there, axis=1),
                                np.sum(tangent_here * tangent_there, axis=1)))
    pitch = np.degrees(np.arcsin(np.clip(-along, -1.0, 1.0)))
    roll = np.degrees(np.arctan2(np.sum(normal_u_here * binormal_there, axis=1),
                                 np.sum(normal_u_here * normal_u_there, axis=1)))
    ok = ok_here & ok_there & ok_normal_here & ok_normal_there
    return np.where(ok, yaw, 0.0), np.where(ok, pitch, 0.0), np.where(ok, roll, 0.0)


def _dihedral(normal_a: np.ndarray, normal_b: np.ndarray, edge: np.ndarray) -> float:
    """Signed fold from ``normal_a`` to ``normal_b`` about ``edge``, in degrees."""
    fold = _dihedral_block(
        np.asarray(normal_a, float).reshape(1, 3), np.asarray(normal_b, float).reshape(1, 3),
        np.asarray(edge, float).reshape(1, 3),
    )
    return float(fold[0])


def _dihedral_block(normal_a: np.ndarray, normal_b: np.ndarray, edge: np.ndarray) -> np.ndarray:
    """``_dihedral`` for one row per edge. A missing direction or normal stays at zero."""
    direction, ok_edge = _unit_rows(edge)
    ok = ok_edge & (np.linalg.norm(normal_a, axis=1) > 1e-15) & (np.linalg.norm(normal_b, axis=1) > 1e-15)
    signed = np.sum(np.cross(normal_a, normal_b) * direction, axis=1)
    fold = np.degrees(np.arctan2(signed, np.sum(normal_a * normal_b, axis=1)))
    return np.where(ok, fold, 0.0)


def _other_face(element_id: int, users: list[int], surfaces) -> int | None:
    """The face that continues this one across an edge, preferring the same surface."""
    others = [i for i in users if i != element_id]
    if not others:
        return None
    same = [i for i in others if surfaces[i] == surfaces[element_id]]
    return same[0] if same else others[0]


def _dihedral_pair(users: list[int], surfaces) -> tuple[int, int] | None:
    """Two faces whose fold is the edge dihedral, lower element id first.

    A same-surface pair wins, so a web landing on a skin edge does not hide the
    fold of the skin. The sign is then from the lower element id to the higher.
    """
    if len(users) < 2:
        return None
    same = [
        (left, right)
        for index, left in enumerate(users)
        for right in users[index + 1:]
        if surfaces[left] == surfaces[right]
    ]
    left, right = same[0] if same else (users[0], users[1])
    return (left, right) if left < right else (right, left)


def _newell_units(points: np.ndarray) -> np.ndarray:
    """Unit Newell normals of equal-sided loops. ``points`` is ``(n, k, 3)``."""
    nxt = np.roll(points, -1, axis=1)
    raw = np.stack([
        np.sum((points[..., 1] - nxt[..., 1]) * (points[..., 2] + nxt[..., 2]), axis=1),
        np.sum((points[..., 2] - nxt[..., 2]) * (points[..., 0] + nxt[..., 0]), axis=1),
        np.sum((points[..., 0] - nxt[..., 0]) * (points[..., 1] + nxt[..., 1]), axis=1),
    ], axis=1)
    length = np.linalg.norm(raw, axis=1)
    unit = np.divide(raw, length[:, None], out=np.zeros_like(raw), where=(length > 1e-15)[:, None])
    return np.where((length > 1e-15)[:, None], unit, raw)


def _corner_angles_block(points: np.ndarray) -> np.ndarray:
    """Corner angles in degrees for equal-sized loops. ``points`` is ``(n, k, 3)``."""
    prev = np.roll(points, 1, axis=1) - points
    nxt = np.roll(points, -1, axis=1) - points
    la = np.linalg.norm(prev, axis=2)
    lb = np.linalg.norm(nxt, axis=2)
    cos = np.sum(prev * nxt, axis=2) / np.maximum(la * lb, 1e-30)
    angles = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    return np.where((la <= 1e-15) | (lb <= 1e-15), 0.0, angles)


def _by_size(elements) -> dict[int, list[int]]:
    grouped: dict[int, list[int]] = {}
    for i, element in enumerate(elements):
        grouped.setdefault(len(element), []).append(i)
    return grouped


def _element_normals(nodes: np.ndarray, elements) -> np.ndarray:
    nodes = np.asarray(nodes, float)
    normals = np.zeros((len(elements), 3), float)
    for n, idx in _by_size(elements).items():
        if n < 3:
            continue
        if n not in (3, 4):
            for i in idx:
                normals[i] = _face_normal(nodes[list(elements[i])])
            continue
        corners = np.asarray([list(elements[i]) for i in idx], np.int64)
        normals[np.asarray(idx, np.int64)] = _newell_units(nodes[corners])
    return normals


def _diagonal_min(nodes: np.ndarray, element) -> float:
    """Shortest distance between nodes of one element that do not share an edge."""
    ids = [int(k) for k in element]
    n = len(ids)
    best = math.inf
    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            gap = float(np.linalg.norm(nodes[ids[i]] - nodes[ids[j]]))
            if gap < best:
                best = gap
    return 0.0 if best is math.inf else best


def diagonal_mins(mesh) -> np.ndarray:
    """One ``diagonal_min`` per element. A quad has two diagonals; a polygon uses the same rule."""
    nodes = np.asarray(mesh.nodes, float)
    out = np.zeros(len(mesh.elements), float)
    for n, idx in _by_size(mesh.elements).items():
        if n < 4:
            continue
        if n != 4:
            for i in idx:
                out[i] = _diagonal_min(nodes, mesh.elements[i])
            continue
        corners = np.asarray([list(mesh.elements[i]) for i in idx], np.int64)
        block = nodes[corners]
        gaps = np.minimum(
            np.linalg.norm(block[:, 0] - block[:, 2], axis=1),
            np.linalg.norm(block[:, 1] - block[:, 3], axis=1),
        )
        out[np.asarray(idx, np.int64)] = gaps
    return out


def _chain_turns(edges: np.ndarray, nodes: np.ndarray, names, normals) -> None:
    """Yaw and pitch from the previous spanwise edge of the same region pair.

    The chain follows increasing ``z``. An edge with no predecessor stays at zero.
    """
    edges["chain_yaw_deg"] = 0.0
    edges["chain_pitch_deg"] = 0.0
    span = np.flatnonzero(edges["kind"] == "spanwise")
    if len(span) == 0:
        return
    name_of = np.asarray([str(n) for n in names], dtype=object)
    groups: dict[tuple[str, ...], list[int]] = {}
    for span_index in span:
        users = _edge_users(edges[span_index])
        pair = tuple(sorted({str(name_of[i]) for i in users}))
        groups.setdefault(pair, []).append(int(span_index))
    query_index: list[int] = []
    prev_rows: list[np.ndarray] = []
    vector_rows: list[np.ndarray] = []
    normal_rows: list[np.ndarray] = []
    for group in groups.values():
        a = edges["a"][group]
        b = edges["b"][group]
        za = nodes[a, 2]
        zb = nodes[b, 2]
        flip = zb + 1e-12 < za
        start = np.where(flip, b, a)
        end = np.where(flip, a, b)
        vectors = np.column_stack([edges["dx"][group], edges["dy"][group], edges["dz"][group]])
        vectors = np.where(flip[:, None], -vectors, vectors)
        order = np.argsort(nodes[start, 2])
        arriving: dict[int, tuple[float, np.ndarray, int]] = {}
        for local in order:
            index = group[int(local)]
            s = int(start[local])
            e = int(end[local])
            vector = vectors[local]
            prev = arriving.get(s)
            if prev is not None:
                _z, prev_vec, prev_index = prev
                users = _edge_users(edges[index])
                normal = normals[users[0]] if users else np.zeros(3)
                if float(np.linalg.norm(normal)) <= 1e-15:
                    prev_users = _edge_users(edges[prev_index])
                    if prev_users:
                        normal = normals[prev_users[0]]
                query_index.append(index)
                prev_rows.append(prev_vec)
                vector_rows.append(vector)
                normal_rows.append(np.asarray(normal, float))
            if e not in arriving:
                arriving[e] = (float(nodes[s, 2]), vector, index)
    if not query_index:
        return
    yaw, pitch, _roll = _trajectory_block(
        np.asarray(prev_rows, float), np.asarray(vector_rows, float),
        np.asarray(normal_rows, float), np.asarray(normal_rows, float),
    )
    chosen = np.asarray(query_index, np.int64)
    edges["chain_yaw_deg"][chosen] = yaw
    edges["chain_pitch_deg"][chosen] = pitch


def build_edges(mesh) -> np.ndarray:
    """One row per undirected node pair, from the element loops."""
    nodes = np.asarray(mesh.nodes, float)
    origins = getattr(mesh, "node_origin", ())
    if origins is None or len(origins) == 0:
        origins = ("station",) * len(nodes)
    origins = np.asarray(origins, dtype=object)
    names = getattr(mesh, "name", ())
    if names is None or len(names) == 0:
        names = ("",) * len(mesh.elements)
    normals = _element_normals(nodes, mesh.elements)
    acc: dict[tuple[int, int], list[int]] = {}
    for i, element in enumerate(mesh.elements):
        for a, b in zip(element, element[1:] + element[:1]):
            if a == b:
                continue
            acc.setdefault((a, b) if a < b else (b, a), []).append(i)
    sag = getattr(mesh, "edge_sagitta", None) or {}
    pairs = list(acc.items())
    out = np.empty(len(pairs), dtype=EDGE_DTYPE)
    if not pairs:
        return out
    ends = np.asarray([key for key, _users in pairs], np.int64)
    a = ends[:, 0]
    b = ends[:, 1]
    vector = nodes[b] - nodes[a]
    left = np.full(len(pairs), -1, np.int64)
    right = np.full(len(pairs), -1, np.int64)
    stored = np.full((len(pairs), 3), -1, np.int64)
    use = np.empty(len(pairs), np.int64)
    meet = np.empty(len(pairs), dtype="U16")
    surface = getattr(mesh, "surface", None)
    if surface is None or len(surface) == 0:
        surface = ("",) * len(mesh.elements)
    for index, ((pa, pb), users) in enumerate(pairs):
        use[index] = len(users)
        stored[index, :min(3, len(users))] = users[:3]
        pair = _dihedral_pair(users, surface)
        if pair is not None:
            left[index], right[index] = pair
        meet[index] = _meet({surface[i] for i in users}, len(users))
    station = (origins[a] == "station") & (origins[b] == "station")
    same_z = np.abs(nodes[a, 2] - nodes[b, 2]) <= 1e-9
    kind = np.full(len(pairs), "interior", dtype="U16")
    kind[station & same_z] = "station"
    kind[station & ~same_z] = "spanwise"
    sagitta = np.fromiter(
        (float(sag.get((int(pa), int(pb)), 0.0)) for (pa, pb), _users in pairs),
        dtype=float, count=len(pairs),
    )
    fold = np.zeros(len(pairs), float)
    paired = left >= 0
    if np.any(paired):
        fold[paired] = _dihedral_block(normals[left[paired]], normals[right[paired]], vector[paired])
    out["a"] = a
    out["b"] = b
    out["use"] = use
    out["e0"] = stored[:, 0]
    out["e1"] = stored[:, 1]
    out["e2"] = stored[:, 2]
    out["length"] = np.linalg.norm(vector, axis=1)
    out["dx"] = vector[:, 0]
    out["dy"] = vector[:, 1]
    out["dz"] = vector[:, 2]
    out["dihedral_deg"] = fold
    out["kind"] = kind
    out["meet"] = meet
    out["sagitta"] = sagitta
    _chain_turns(out, nodes, names, normals)
    return out


def _chain_open(edges: np.ndarray, nodes: np.ndarray, z_root: float) -> tuple[OpenLoop, ...]:
    open_edges = edges[edges["use"] == 1] if len(edges) else edges
    unused = {(int(edge["a"]), int(edge["b"])) for edge in open_edges}
    by_node: dict[int, list[tuple[int, int]]] = {}
    for a, b in unused:
        by_node.setdefault(a, []).append((a, b))
        by_node.setdefault(b, []).append((a, b))
    loops: list[OpenLoop] = []
    seen: set[tuple[int, int]] = set()
    for start in list(unused):
        if start in seen:
            continue
        ring = [start[0]]
        cur = start[1]
        seen.add(start)
        ring.append(cur)
        guard = 0
        while cur != ring[0] and guard < len(unused) + 2:
            guard += 1
            nxt = None
            for edge in by_node.get(cur, ()):
                if edge in seen:
                    continue
                nxt = edge[0] if edge[1] == cur else edge[1]
                seen.add(edge)
                break
            if nxt is None:
                break
            if nxt == ring[0]:
                break
            ring.append(nxt)
            cur = nxt
        zs = nodes[ring, 2]
        z = float(np.median(zs))
        is_root = bool(np.all(np.abs(zs - z_root) <= 1e-8))
        loops.append(OpenLoop(z=z, n_edges=max(len(ring) - (0 if cur != ring[0] else 0), 1),
                              is_root=is_root, nodes=tuple(int(k) for k in ring)))
    return tuple(loops)


def valence_per_node(mesh, edges: np.ndarray) -> np.ndarray:
    """Skin-edge valence when the node is on the skin, otherwise the valence on its own surface."""
    n_nodes = len(mesh.nodes)
    skin_node = np.zeros(n_nodes, bool)
    surfaces = np.asarray(list(mesh.surface), dtype=object)
    for element, surf in zip(mesh.elements, surfaces):
        if surf == "skin":
            for node in element:
                skin_node[int(node)] = True
    count = np.zeros(n_nodes, np.int64)
    if len(edges) == 0:
        return count
    a = np.asarray(edges["a"], np.int64)
    b = np.asarray(edges["b"], np.int64)
    skin_edge = np.zeros(len(edges), bool)
    for name in ("e0", "e1", "e2"):
        ids = np.asarray(edges[name], np.int64)
        ok = ids >= 0
        if np.any(ok):
            skin_edge[ok] |= surfaces[ids[ok]] == "skin"
    take_a = np.where(skin_node[a], skin_edge, True)
    take_b = np.where(skin_node[b], skin_edge, True)
    np.add.at(count, a, take_a.astype(np.int64))
    np.add.at(count, b, take_b.astype(np.int64))
    return count


def _cell_groups(key: np.ndarray) -> dict[tuple[int, int, int], np.ndarray]:
    """Node indices in each integer cell. ``key`` is ``(n, 3)``."""
    if len(key) == 0:
        return {}
    order = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    sorted_key = key[order]
    change = np.any(sorted_key[1:] != sorted_key[:-1], axis=1) if len(order) > 1 else np.zeros(0, dtype=bool)
    starts = np.concatenate([[0], np.flatnonzero(change) + 1])
    stops = np.concatenate([starts[1:], [len(order)]])
    return {
        (int(sorted_key[start, 0]), int(sorted_key[start, 1]), int(sorted_key[start, 2])): order[start:stop]
        for start, stop in zip(starts, stops)
    }


def near_misses(nodes: np.ndarray, tol: float = 1e-8) -> tuple[Seam, ...]:
    """Distinct node ids closer than ``tol``. The mesh index merges exact coordinates only."""
    pts = np.asarray(nodes, float)
    if len(pts) == 0:
        return ()
    cell = max(tol * 2.0, 1e-9)
    groups = _cell_groups(np.floor(pts / cell).astype(np.int64))
    left: list[np.ndarray] = []
    right: list[np.ndarray] = []
    for base, idx in groups.items():
        for dz in (0, 1):
            for dy in (0, 1):
                for dx in (0, 1):
                    other = groups.get((base[0] + dx, base[1] + dy, base[2] + dz))
                    if other is None:
                        continue
                    if dx == dy == dz == 0:
                        if len(idx) < 2:
                            continue
                        ii, jj = np.triu_indices(len(idx), k=1)
                        a = idx[ii]
                        b = idx[jj]
                    else:
                        a = np.repeat(idx, len(other))
                        b = np.tile(other, len(idx))
                    keep = b > a
                    if np.any(keep):
                        left.append(a[keep])
                        right.append(b[keep])
    if not left:
        return ()
    ia = np.concatenate(left)
    ib = np.concatenate(right)
    dist = np.linalg.norm(pts[ib] - pts[ia], axis=1)
    hit = (dist > 0.0) & (dist < tol)
    return tuple(
        Seam("near_miss", int(i), int(j), int(j), float(d))
        for i, j, d in zip(ia[hit], ib[hit], dist[hit])
    )


def _segment_hits(pts: np.ndarray, edge_ids: np.ndarray, node_ids: np.ndarray,
                  a: np.ndarray, b: np.ndarray, pa: np.ndarray, ab: np.ndarray,
                  span2: np.ndarray, lo: np.ndarray, hi: np.ndarray, tol: float
                  ) -> list[tuple[int, int, int, float]]:
    """Nodes of ``node_ids`` that lie in the interior of ``edge_ids``."""
    if len(edge_ids) == 0 or len(node_ids) == 0:
        return []
    found: list[tuple[int, int, int, float]] = []
    step = max(1, 200000 // max(len(node_ids), 1))
    for start in range(0, len(edge_ids), step):
        batch = edge_ids[start:start + step]
        ee = np.repeat(batch, len(node_ids))
        nn = np.tile(node_ids, len(batch))
        point = pts[nn]
        ok = np.all(point >= lo[ee], axis=1) & np.all(point <= hi[ee], axis=1)
        ok &= (nn != a[ee]) & (nn != b[ee])
        if not np.any(ok):
            continue
        ee = ee[ok]
        nn = nn[ok]
        point = point[ok]
        t = np.sum((point - pa[ee]) * ab[ee], axis=1) / span2[ee]
        interior = (t > 1e-4) & (t < 1.0 - 1e-4)
        if not np.any(interior):
            continue
        ee = ee[interior]
        nn = nn[interior]
        point = point[interior]
        t = t[interior]
        dist = np.linalg.norm(point - (pa[ee] + t[:, None] * ab[ee]), axis=1)
        close = dist <= tol
        for edge, node, gap in zip(ee[close], nn[close], dist[close]):
            found.append((int(a[edge]), int(b[edge]), int(node), float(gap)))
    return found


def t_junctions(nodes: np.ndarray, edges: np.ndarray, tol: float = 1e-5) -> tuple[Seam, ...]:
    """A node lying in an edge without being one of its ends."""
    pts = np.asarray(nodes, float)
    if len(edges) == 0 or len(pts) == 0:
        return ()
    cell = 0.05
    length = np.asarray(edges["length"], float)
    live = np.flatnonzero(length > tol)
    if len(live) == 0:
        return ()
    a = np.asarray(edges["a"], np.int64)[live]
    b = np.asarray(edges["b"], np.int64)[live]
    pa = pts[a]
    pb = pts[b]
    ab = pb - pa
    span2 = np.maximum(length[live] ** 2, 1e-30)
    lo = np.minimum(pa, pb) - tol
    hi = np.maximum(pa, pb) + tol
    groups = _cell_groups(np.floor(pts / cell).astype(np.int64))
    klo = np.floor(lo / cell).astype(np.int64)
    khi = np.floor(hi / cell).astype(np.int64)
    by_cell: dict[tuple[int, int, int], list[int]] = {}
    for edge in range(len(live)):
        for ix in range(int(klo[edge, 0]), int(khi[edge, 0]) + 1):
            for iy in range(int(klo[edge, 1]), int(khi[edge, 1]) + 1):
                for iz in range(int(klo[edge, 2]), int(khi[edge, 2]) + 1):
                    by_cell.setdefault((ix, iy, iz), []).append(edge)
    found: list[tuple[int, int, int, float]] = []
    for key, edge_ids in by_cell.items():
        nodes_here = groups.get(key)
        if nodes_here is None:
            continue
        found.extend(_segment_hits(
            pts, np.asarray(edge_ids, np.int64), nodes_here, a, b, pa, ab, span2, lo, hi, tol,
        ))
    return tuple(Seam("t_junction", ia, ib, node, dist) for ia, ib, node, dist in found)


def _ruled_valence_defects(mesh, edges: np.ndarray, valence: Sequence[int] | np.ndarray) -> tuple[int, ...]:
    """Interior skin nodes of a ruled bay whose valence is not 4.

    A skin edge used by one skin element is the boundary of the sheet, including
    where a rib or the tip closure covers it. Those nodes are not interior.
    """
    skin_use: dict[tuple[int, int], int] = {}
    for element, surface in zip(mesh.elements, mesh.surface):
        if surface != "skin":
            continue
        for a, b in zip(element, element[1:] + element[:1]):
            key = (a, b) if a < b else (b, a)
            skin_use[key] = skin_use.get(key, 0) + 1
    on_boundary = {n for key, count in skin_use.items() if count == 1 for n in key}
    ruled_nodes: dict[int, bool] = {}
    for element, surface, constructor in zip(mesh.elements, mesh.surface, mesh.constructor):
        if surface != "skin":
            continue
        for node in element:
            if constructor != "ruled":
                ruled_nodes[node] = False
            else:
                ruled_nodes.setdefault(node, True)
    return tuple(
        node for node, only_ruled in ruled_nodes.items()
        if only_ruled and node not in on_boundary and valence[node] not in (0, 4)
    )


def trailing_face_failures(mesh, *, angle: float = 20.0, diagonal: float = 1e-3,
                          edge: float = 5e-3) -> list[dict]:
    """Trailing-wrap quads that are shallow, inverted, or collapsed across a real edge.

    Adapter over :func:`~implicit_geometry_kernel.mesh.health.run_mesh_health`.
    A row is region ``te``. Reasons are ``angle`` / ``inverted`` / ``diagonal``.
    """
    from implicit_geometry_kernel.mesh.health import run_mesh_health

    reason_of = {
        "trailing_face_angle": "angle",
        "trailing_face_inverted": "inverted",
        "trailing_face_diagonal": "diagonal",
    }
    by_element: dict[int, list[str]] = {}
    meta: dict[int, tuple[int, str]] = {}
    for hit in run_mesh_health(
            mesh, angle=angle, trailing_diagonal=diagonal, trailing_edge=edge,
            include_geometry=False):
        reason = reason_of.get(hit.code)
        if reason is None or hit.element < 0:
            continue
        by_element.setdefault(hit.element, []).append(reason)
        meta[hit.element] = (hit.bay, hit.constructor)
    rows = []
    for element, reasons in by_element.items():
        bay, constructor = meta[element]
        rows.append({
            "element": element,
            "bay": bay,
            "constructor": constructor,
            "reasons": tuple(dict.fromkeys(reasons)),
        })
    return rows


def transition_band_failures(mesh, *, angle: float = 20.0, area_ratio: float = 2.0) -> list[dict]:
    """Skin and web quads that kink, or jump in area, inside one region and bay.

    Adapter over :func:`~implicit_geometry_kernel.mesh.health.run_mesh_health`.
    """
    from implicit_geometry_kernel.mesh.health import run_mesh_health

    reason_of = {
        "transition_band_angle": "angle",
        "transition_band_area": "area",
    }
    by_element: dict[int, list[str]] = {}
    meta: dict[int, tuple[str, int, str]] = {}
    for hit in run_mesh_health(
            mesh, angle=angle, area_ratio=area_ratio, include_geometry=False):
        reason = reason_of.get(hit.code)
        if reason is None or hit.element < 0:
            continue
        by_element.setdefault(hit.element, []).append(reason)
        meta[hit.element] = (hit.region, hit.bay, hit.constructor)
    rows = []
    for element, reasons in by_element.items():
        region, bay, constructor = meta[element]
        rows.append({
            "element": element,
            "region": region,
            "bay": bay,
            "constructor": constructor,
            "reasons": tuple(dict.fromkeys(reasons)),
        })
    return rows


def mesh_failures(mesh, *, angle: float = 20.0, skew: float | None = None,
                  volume_ratio: float | None = None, sagitta: float | None = None,
                  distance: float | None = None) -> list[dict]:
    """Elements past a threshold, grouped by constructor and bay. Not one count.

    Adapter over :func:`~implicit_geometry_kernel.mesh.health.run_mesh_health`.
    """
    from implicit_geometry_kernel.mesh.health import run_mesh_health

    quality_codes = {
        "quad_angle", "quad_skew", "quad_volume_ratio", "quad_inverted",
        "section_distance", "edge_sagitta",
    }
    hits_by_key: dict[tuple, list[int]] = {}
    for hit in run_mesh_health(
            mesh, angle=angle, skew=skew, volume_ratio=volume_ratio,
            sagitta=sagitta, distance=distance, include_geometry=False):
        if hit.code not in quality_codes or hit.element < 0:
            continue
        key = (hit.constructor, hit.bay)
        hits_by_key.setdefault(key, []).append(hit.element)
    worst_at: dict[int, np.void] = {}
    corner_rows = getattr(mesh, "corners", None)
    if corner_rows is not None and getattr(getattr(corner_rows, "dtype", None), "names", None):
        for corner in corner_rows:
            element_id = int(corner["element_id"])
            kept = worst_at.get(element_id)
            if kept is None or float(corner["angle_deg"]) < float(kept["angle_deg"]):
                worst_at[element_id] = corner
    groups = []
    for (constructor, bay), hit_elements in hits_by_key.items():
        elements = tuple(dict.fromkeys(hit_elements))
        z0 = float(mesh.z0[elements[0]]) if mesh.z0 else 0.0
        z1 = float(mesh.z1[elements[0]]) if mesh.z1 else 0.0
        worst = None
        for element_id in elements:
            corner = worst_at.get(element_id)
            if corner is not None and (worst is None or float(corner["angle_deg"]) < float(worst["angle_deg"])):
                worst = corner
        groups.append({
            "constructor": constructor, "bay": bay, "z0": z0, "z1": z1,
            "elements": elements, "n": len(elements),
            "corner": None if worst is None else int(worst["corner"]),
        })
    return groups


def _csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _csv_columns(path: Path, header: list[str], columns: list) -> None:
    """Write one CSV from parallel columns. An empty table is an empty file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(columns[0]) if columns else 0
    if n == 0:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(zip(*columns))


def off_root_open_edges(mesh, edges: np.ndarray | None = None) -> tuple[dict, ...]:
    """Open edges whose nodes are not both on the root station, grouped by region and kind.

    Each group is ``{"region", "kind", "count", "z0", "z1"}`` sorted by count
    descending. Empty when the only free edges lie on the root ring.
    """
    from implicit_geometry_kernel.mesh.health import run_mesh_health

    # Prefer registry hits when the mesh already carries assembled edges.
    if edges is None and getattr(mesh, "edges", None) is not None and len(mesh.edges):
        groups: dict[tuple[str, str], dict] = {}
        edge_rows = np.asarray(mesh.edges)
        for hit in run_mesh_health(mesh):
            if hit.code != "open_edge_off_root" or hit.edge < 0:
                continue
            kind = str(edge_rows["kind"][hit.edge]) if hit.edge < len(edge_rows) else hit.detail
            key = (hit.region, kind)
            row = groups.get(key)
            if row is None:
                groups[key] = {
                    "region": hit.region, "kind": kind, "count": 1,
                    "z0": float(hit.z0), "z1": float(hit.z1),
                }
            else:
                row["count"] += 1
                row["z0"] = min(row["z0"], float(hit.z0))
                row["z1"] = max(row["z1"], float(hit.z1))
        return tuple(sorted(
            groups.values(),
            key=lambda row: (-int(row["count"]), row["region"], row["kind"]),
        ))
    edges = build_edges(mesh) if edges is None else np.asarray(edges)
    nodes = np.asarray(mesh.nodes, float)
    if len(edges) == 0 or len(nodes) == 0:
        return ()
    root_z = float(np.min(nodes[:, 2]))
    names = mesh.name if len(mesh.name) == len(mesh.elements) else ("",) * len(mesh.elements)
    groups = {}
    for edge in edges:
        if int(edge["use"]) != 1:
            continue
        a, b = int(edge["a"]), int(edge["b"])
        za, zb = float(nodes[a, 2]), float(nodes[b, 2])
        if abs(za - root_z) <= 1e-9 and abs(zb - root_z) <= 1e-9:
            continue
        el = int(edge["e0"])
        region = str(names[el]) if 0 <= el < len(names) else ""
        kind = str(edge["kind"])
        key = (region, kind)
        row = groups.get(key)
        z_lo, z_hi = min(za, zb), max(za, zb)
        if row is None:
            groups[key] = {"region": region, "kind": kind, "count": 1, "z0": z_lo, "z1": z_hi}
        else:
            row["count"] += 1
            row["z0"] = min(row["z0"], z_lo)
            row["z1"] = max(row["z1"], z_hi)
    ordered = sorted(groups.values(), key=lambda row: (-int(row["count"]), row["region"], row["kind"]))
    return tuple(ordered)


def _summary(mesh, hits=None) -> str:
    from implicit_geometry_kernel.mesh.health import format_health_summary, run_mesh_health

    table = np.asarray(mesh.element_quality)
    quads = [i for i, e in enumerate(mesh.elements) if len(e) == 4]
    under = 0
    inverted = []
    if len(table) == len(mesh.elements) and quads:
        under = int(np.sum(table["angle_min"][quads] < 20.0))
        inverted = [i for i in range(len(mesh.elements)) if float(table["jacobian_min"][i]) < 0.0]
    ratio = float(np.max(table["volume_ratio"])) if len(table) else 1.0
    leaks = [lp for lp in mesh.loops if not lp.is_root]
    if hits is None:
        hits = run_mesh_health(mesh)
    off_n = sum(1 for hit in hits if hit.code == "open_edge_off_root")
    lines = [
        f"quads under 20 deg: {under}  (elements.csv angle_min)",
        f"volume ratio max: {ratio:.4f}  (elements.csv volume_ratio)",
        f"open loops that are not the root: {len(leaks)}  (loops.csv)",
        f"open edges off the root: {off_n}",
        f"t-junctions: {len(mesh.junctions)}  (seams are junctions.csv)",
        f"inverted elements: {len(inverted)}  (elements.csv jacobian_min)",
        f"near misses: {len(mesh.near_miss)}",
        f"ruled valence defects: {len(mesh.valence_defects)}",
    ]
    lines.extend(format_health_summary(hits))
    if leaks:
        lines.append("leak rows: " + ", ".join(str(i) for i, lp in enumerate(mesh.loops) if not lp.is_root))
    if inverted:
        lines.append("inverted element ids: " + ", ".join(str(i) for i in inverted[:40]))
    shallow = _shallowest_quad_corner(mesh)
    if shallow is not None:
        lines.append(
            f"shallowest quad corner: element {int(shallow['element_id'])} corner {int(shallow['corner'])} "
            f"node {int(shallow['node_id'])} xyz {float(shallow['x']):.6g} {float(shallow['y']):.6g} "
            f"{float(shallow['z']):.6g} "
            f"angle {float(shallow['angle_deg']):.4f} deg constructor {shallow['constructor']} "
            f"bay {int(shallow['bay'])} region {shallow['region']}"
        )
    return "\n".join(lines) + "\n"


def _shallowest_quad_corner(mesh):
    """The quad corner with the smallest interior angle, or ``None`` when the mesh has no quad corner."""
    best = None
    elements = mesh.elements
    corners = getattr(mesh, "corners", None)
    if corners is None or not getattr(getattr(corners, "dtype", None), "names", None):
        return None
    for corner in corners:
        element_id = int(corner["element_id"])
        if element_id >= len(elements) or len(elements[element_id]) != 4:
            continue
        if best is None or float(corner["angle_deg"]) < float(best["angle_deg"]):
            best = corner
    return best


def build_corners(mesh, edges: np.ndarray | None = None) -> np.ndarray:
    """One row per corner of every element. The angle is the one :func:`_corner_angles` returns at that loop index."""
    nodes = np.asarray(mesh.nodes, float)
    n_elements = len(mesh.elements)
    names = mesh.name if len(mesh.name) == n_elements else ("",) * n_elements
    surfaces = mesh.surface if len(mesh.surface) == n_elements else ("",) * n_elements
    constructors = mesh.constructor if len(mesh.constructor) == n_elements else ("",) * n_elements
    bays = mesh.bay if len(mesh.bay) == n_elements else (-1,) * n_elements
    normals = _element_normals(nodes, mesh.elements)
    angle_of: list[np.ndarray | None] = [None] * n_elements
    for n, idx in _by_size(mesh.elements).items():
        if n not in (3, 4):
            continue
        corners = np.asarray([list(mesh.elements[i]) for i in idx], np.int64)
        block = _corner_angles_block(nodes[corners])
        for row, i in zip(block, idx):
            angle_of[i] = row
    if edges is None:
        users: dict[tuple[int, int], list[int]] = {}
        for element_id, element in enumerate(mesh.elements):
            ids = [int(k) for k in element]
            for a, b in zip(ids, ids[1:] + ids[:1]):
                if a == b:
                    continue
                users.setdefault((a, b) if a < b else (b, a), []).append(element_id)
    else:
        users = {
            (int(edge["a"]), int(edge["b"])): _edge_users(edge)
            for edge in edges
        }
    n_rows = sum(len(element) for element in mesh.elements if len(element) >= 3)
    rows = np.empty(n_rows, dtype=CORNER_DTYPE)
    if n_rows == 0:
        return rows
    element_ids = np.empty(n_rows, np.int64)
    corner_index = np.empty(n_rows, np.int64)
    node_ids = np.empty(n_rows, np.int64)
    prev_ids = np.empty(n_rows, np.int64)
    next_ids = np.empty(n_rows, np.int64)
    angle = np.empty(n_rows, float)
    neighbor_ids = np.full(n_rows, -1, np.int64)
    cursor = 0
    for element_id, element in enumerate(mesh.elements):
        ids = [int(k) for k in element]
        n = len(ids)
        if n < 3:
            continue
        angles = angle_of[element_id]
        if angles is None:
            angles = _corner_angles(nodes[ids])
        stop = cursor + n
        element_ids[cursor:stop] = element_id
        corner_index[cursor:stop] = np.arange(n)
        node_ids[cursor:stop] = ids
        prev_ids[cursor:stop] = [ids[-1], *ids[:-1]]
        next_ids[cursor:stop] = [*ids[1:], ids[0]]
        angle[cursor:stop] = angles
        for k, node_id in enumerate(ids):
            nxt = next_ids[cursor + k]
            key = (node_id, int(nxt)) if node_id < nxt else (int(nxt), node_id)
            neighbor = _other_face(element_id, users.get(key, []), surfaces)
            if neighbor is not None:
                neighbor_ids[cursor + k] = neighbor
        cursor = stop
    point = nodes[node_ids]
    incoming = point - nodes[prev_ids]
    outgoing = nodes[next_ids] - point
    normal_here = normals[element_ids]
    normal_next = normal_here.copy()
    shared = neighbor_ids >= 0
    normal_next[shared] = normals[neighbor_ids[shared]]
    yaw, pitch, roll = _trajectory_block(incoming, outgoing, normal_here, normal_next)
    rows["element_id"] = element_ids
    rows["corner"] = corner_index
    rows["node_id"] = node_ids
    rows["x"] = point[:, 0]
    rows["y"] = point[:, 1]
    rows["z"] = point[:, 2]
    rows["angle_deg"] = angle
    rows["yaw_deg"] = yaw
    rows["pitch_deg"] = pitch
    rows["roll_deg"] = roll
    rows["dx_in"] = incoming[:, 0]
    rows["dy_in"] = incoming[:, 1]
    rows["dz_in"] = incoming[:, 2]
    rows["dx_out"] = outgoing[:, 0]
    rows["dy_out"] = outgoing[:, 1]
    rows["dz_out"] = outgoing[:, 2]
    rows["prev_id"] = prev_ids
    rows["next_id"] = next_ids
    rows["constructor"] = np.asarray(constructors, dtype="U32")[element_ids]
    rows["bay"] = np.asarray(bays, np.int64)[element_ids]
    rows["region"] = np.asarray(names, dtype="U64")[element_ids]
    rows["surface"] = np.asarray(surfaces, dtype="U16")[element_ids]
    return rows


def assemble(mesh, z_root: float) -> dict:
    """Edges, loops, corners, valence, near-misses, and T-junctions for a mesh whose elements are already set."""
    edges = build_edges(mesh)
    nodes = np.asarray(mesh.nodes, float)
    valence = valence_per_node(mesh, edges)
    return {
        "edges": edges,
        "loops": _chain_open(edges, nodes, float(z_root)),
        "corners": build_corners(mesh, edges),
        "valence": valence,
        "near_miss": near_misses(nodes),
        "junctions": t_junctions(nodes, edges),
        "valence_defects": _ruled_valence_defects(mesh, edges, valence),
        "diagonal_min": diagonal_mins(mesh),
    }


def write_mesh_tables(mesh, directory) -> Path:
    """Write one CSV per mesh table, plus ``summary.txt``. Returns the directory."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    nodes = np.asarray(mesh.nodes, float)
    dist = np.asarray(mesh.section_distance, float) if len(mesh.section_distance) else np.zeros(len(nodes))
    valence = mesh.valence if len(mesh.valence) else (0,) * len(nodes)
    origin = mesh.node_origin if len(mesh.node_origin) else ("station",) * len(nodes)
    node_header = ["node_id", "x", "y", "z", "origin", "section_distance", "source", "valence"]
    node_columns = [
        np.arange(len(nodes)), nodes[:, 0], nodes[:, 1], nodes[:, 2],
        list(origin), dist, (np.arange(len(nodes)) < int(mesh.n_source)).astype(np.int64),
        np.asarray(valence, np.int64),
    ]
    quality = np.asarray(mesh.element_quality)
    n_elements = len(mesh.elements)
    width = max((len(e) for e in mesh.elements), default=4)
    element_header = ["element_id", *[f"n{k}" for k in range(width)],
                      "surface", "region", "partition", "bay", "z0", "z1", "count0", "count1",
                      "constructor", "flipped", "diagonal_min"]
    if len(quality) == n_elements and quality.dtype.names:
        element_header.extend(quality.dtype.names)
    node_cols: list[list[int | str]] = [[""] * n_elements for _ in range(width)]
    for i, element in enumerate(mesh.elements):
        for k, node in enumerate(element):
            node_cols[k][i] = int(node)
    has_bay = len(mesh.bay)
    has_diag = len(mesh.diagonal_min)
    element_columns: list = [
        np.arange(n_elements),
        *node_cols,
        list(mesh.surface) if n_elements else [],
        list(mesh.name) if n_elements else [],
        list(mesh.partition) if n_elements else [],
        np.asarray(mesh.bay, np.int64) if has_bay else np.full(n_elements, -1, np.int64),
        np.asarray(mesh.z0, float) if len(mesh.z0) else np.zeros(n_elements),
        np.asarray(mesh.z1, float) if len(mesh.z1) else np.zeros(n_elements),
        np.asarray(mesh.n0, np.int64) if len(mesh.n0) else np.zeros(n_elements, np.int64),
        np.asarray(mesh.n1, np.int64) if len(mesh.n1) else np.zeros(n_elements, np.int64),
        list(mesh.constructor) if len(mesh.constructor) else [""] * n_elements,
        np.asarray(mesh.flipped, np.int64) if len(mesh.flipped) else np.zeros(n_elements, np.int64),
        np.asarray(mesh.diagonal_min, float) if has_diag else [""] * n_elements,
    ]
    if len(quality) == n_elements and quality.dtype.names:
        element_columns.extend(quality[field] for field in quality.dtype.names)
    edges = np.asarray(mesh.edges)
    edge_header = ["a", "b", "use", "elements", "length", "dx", "dy", "dz", "dihedral_deg",
                   "chain_yaw_deg", "chain_pitch_deg", "kind", "meet", "sagitta"]
    if len(edges):
        e0 = np.asarray(edges["e0"], np.int64)
        e1 = np.asarray(edges["e1"], np.int64)
        e2 = np.asarray(edges["e2"], np.int64)
        user_text = np.empty(len(edges), dtype=object)
        for i in range(len(edges)):
            ids = [e0[i], e1[i], e2[i]]
            user_text[i] = " ".join(str(int(v)) for v in ids if int(v) >= 0)
        edge_columns = [
            edges["a"], edges["b"], edges["use"], user_text,
            edges["length"], edges["dx"], edges["dy"], edges["dz"], edges["dihedral_deg"],
            edges["chain_yaw_deg"], edges["chain_pitch_deg"],
            edges["kind"], edges["meet"], edges["sagitta"],
        ]
    else:
        edge_columns = [np.asarray([]) for _ in edge_header]
    count_rows = [{
        "z": c.z, "region": c.region, "length": c.length, "edge": c.edge,
        "edge_from": c.edge_from, "previous": c.previous, "count": c.count, "reason": c.reason,
    } for c in mesh.counts]
    station_rows = [{
        "z": s.z, "origin": s.origin, "refused": int(s.refused), "detail": s.detail,
    } for s in mesh.station_log]
    strip_rows = [{
        "bay": s.bay, "z0": s.z0, "z1": s.z1, "region": s.region, "n0": s.n0, "n1": s.n1,
        "constructor": s.constructor, "trigger": s.trigger, "adv_a": s.adv_a,
        "adv_b": s.adv_b, "ties": s.ties, "note": s.note, "chord_ratio": s.chord_ratio,
        "site0": "" if s.site0 < 0 else s.site0, "site1": "" if s.site1 < 0 else s.site1,
    } for s in mesh.strips]
    corners = np.asarray(mesh.corners)
    corner_header = ["element_id", "corner", "node_id", "x", "y", "z", "angle_deg",
                     "yaw_deg", "pitch_deg", "roll_deg", "dx_in", "dy_in", "dz_in",
                     "dx_out", "dy_out", "dz_out", "prev_id", "next_id",
                     "constructor", "bay", "region", "surface"]
    corner_columns = [corners[name] for name in corner_header] if len(corners) else [np.asarray([]) for _ in corner_header]
    loop_rows = [{
        "z": lp.z, "n_edges": lp.n_edges, "is_root": int(lp.is_root),
        "nodes": " ".join(str(k) for k in lp.nodes),
    } for lp in mesh.loops]
    seam_rows = [{
        "kind": s.kind, "a": s.a, "b": s.b, "node": s.node, "distance": s.distance,
    } for s in (*mesh.near_miss, *mesh.junctions)]
    _csv_columns(directory / "nodes.csv", node_header, node_columns)
    _csv_columns(directory / "elements.csv", element_header, element_columns)
    _csv_columns(directory / "edges.csv", edge_header, edge_columns)
    _csv(directory / "counts.csv", count_rows)
    _csv(directory / "stations.csv", station_rows)
    _csv(directory / "strips.csv", strip_rows)
    _csv_columns(directory / "corners.csv", corner_header, corner_columns)
    _csv(directory / "loops.csv", loop_rows)
    _csv(directory / "junctions.csv", seam_rows)
    from implicit_geometry_kernel.mesh.health import (
        health_counts_table, hits_to_rows, run_mesh_health,
    )
    hits = run_mesh_health(mesh)
    _csv(directory / "failures.csv", hits_to_rows(hits))
    _csv(directory / "health_counts.csv", health_counts_table(hits))
    (directory / "summary.txt").write_text(_summary(mesh, hits), encoding="utf-8")
    return directory


__all__ = [
    "CONSTRUCTORS", "CORNER_DTYPE", "EDGE_DTYPE", "KINDS", "MEETS", "ORIGINS",
    "Corner", "MeshEdge", "OpenLoop", "Seam", "StationDecision", "StripDecision",
    "assemble", "build_corners", "build_edges", "mesh_failures", "near_misses",
    "off_root_open_edges", "t_junctions", "trailing_face_failures",
    "transition_band_failures", "valence_per_node", "write_mesh_tables",
]
