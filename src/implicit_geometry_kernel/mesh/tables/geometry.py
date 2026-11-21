"""Edge and corner tables: dihedrals, chain turns, element normals and corner angles.

Part of :mod:`implicit_geometry_kernel.mesh.tables.tables`, split out of it in 0.9.4.
"""
from __future__ import annotations

import math

import numpy as np

from implicit_geometry_kernel.mesh.quality import _corner_angles
from implicit_geometry_kernel.mesh.tables.records import (
    CORNER_DTYPE,
    EDGE_DTYPE,
    _edge_users,
)


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
