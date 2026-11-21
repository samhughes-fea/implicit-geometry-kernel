"""CSV dump of the mesh tables and the summary written with it.

Part of :mod:`implicit_geometry_kernel.mesh.tables.tables`, split out of it in 0.9.4.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from implicit_geometry_kernel.mesh.tables.geometry import (
    build_edges,
)


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
