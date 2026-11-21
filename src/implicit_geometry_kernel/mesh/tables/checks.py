"""Mesh failure queries: open loops, valence, near misses, T-junctions and the face checks.

Part of :mod:`implicit_geometry_kernel.mesh.tables.tables`, split out of it in 0.9.4.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from implicit_geometry_kernel.mesh.tables.records import (
    OpenLoop,
    Seam,
)


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
