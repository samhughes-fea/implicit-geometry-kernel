"""Write a :class:`~implicit_geometry_kernel.mesh.shell.ShellMesh` as legacy ASCII VTK.

VTK writes each cell as stored: triangles (``VTK_TRIANGLE``), quads
(``VTK_QUAD``), and any remaining polygon (``VTK_POLYGON``). Cell data: ``thickness``,
``eccentricity``, the element-quality columns, and integer ids for ``surface``
and ``partition`` (their names are listed in the file header comment order,
see ``partition_names``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from implicit_geometry_kernel.mesh.health import (
    HEALTH_CODES, HEALTH_SEVERITIES, edge_health_codes, element_health_codes,
    run_mesh_health,
)
from implicit_geometry_kernel.mesh.shell import SURFACES, ShellMesh
from implicit_geometry_kernel.mesh.tables import CONSTRUCTORS, KINDS, MEETS, ORIGINS

_VTK_VERTEX, _VTK_LINE, _VTK_TRIANGLE, _VTK_QUAD, _VTK_POLYGON = 1, 3, 5, 9, 7


def _vtk_cell(element) -> int:
    n = len(element)
    if n == 3:
        return _VTK_TRIANGLE
    if n == 4:
        return _VTK_QUAD
    return _VTK_POLYGON


def _scalars(handle, name: str, kind: str, values) -> None:
    handle.write(f"SCALARS {name} {kind} 1\nLOOKUP_TABLE default\n")
    if kind == "int":
        np.savetxt(handle, np.asarray(values, np.int64).reshape(-1), fmt="%d")
    else:
        np.savetxt(handle, np.asarray(values, float).reshape(-1), fmt="%.16g")


def _codes(names: tuple[str, ...], values) -> np.ndarray:
    lookup = {name: i for i, name in enumerate(names)}
    fallback = len(names)
    return np.fromiter((lookup.get(str(value), fallback) for value in values), np.int64, count=len(values))


def _legend(mesh: ShellMesh, partitions: tuple[str, ...]) -> str:
    regions = tuple(dict.fromkeys(mesh.name))
    parts = [
        "surface " + " ".join(f"{i}={name}" for i, name in enumerate(SURFACES)),
        "origin " + " ".join(f"{i}={name}" for i, name in enumerate(ORIGINS)),
        "constructor " + " ".join(f"{i}={name}" for i, name in enumerate(CONSTRUCTORS)),
        "kind " + " ".join(f"{i}={name}" for i, name in enumerate(KINDS)),
        "meet " + " ".join(f"{i}={name}" for i, name in enumerate(MEETS)),
        "health_code " + " ".join(f"{i}={name}" for i, name in enumerate(HEALTH_CODES))
        + f" {len(HEALTH_CODES)}=none",
        "health_severity " + " ".join(f"{i}={name}" for i, name in enumerate(HEALTH_SEVERITIES))
        + f" {len(HEALTH_SEVERITIES)}=none",
        "region " + " ".join(f"{i}={name}" for i, name in enumerate(regions)),
        "partition " + " ".join(f"{i}={name}" for i, name in enumerate(partitions)),
    ]
    return "; ".join(parts)


def write_vtk(mesh: ShellMesh, path: str | Path, title: str = "shell mesh") -> tuple[str, ...]:
    """Write ``mesh`` to ``path``; returns the partition names in id order.

    Face cells stay in this file. Edge lines go to ``<stem>.edges.vtk`` through
    the same point indices. The integer legends are packed into the title line
    and written in full to ``<stem>.legend.txt``.
    """
    path = Path(path)
    partitions = tuple(dict.fromkeys(mesh.partition))
    regions = tuple(dict.fromkeys(mesh.name))
    n_cells = len(mesh.elements)
    n_points = len(mesh.nodes)
    size = sum(len(e) + 1 for e in mesh.elements)
    legend = _legend(mesh, partitions)
    head = " ".join(str(title).split()) or "shell mesh"
    nodes = np.asarray(mesh.nodes, float)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write((head + " | " + legend)[:255] + "\n")
        handle.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        handle.write(f"POINTS {n_points} double\n")
        np.savetxt(handle, nodes, fmt="%.16g")
        handle.write(f"CELLS {n_cells} {size}\n")
        sizes = np.fromiter((len(e) for e in mesh.elements), np.int64, count=n_cells)
        max_n = int(sizes.max()) if n_cells else 0
        if n_cells and np.all(sizes == sizes[0]):
            block = np.empty((n_cells, 1 + int(sizes[0])), np.int64)
            block[:, 0] = sizes
            block[:, 1:] = np.asarray([list(e) for e in mesh.elements], np.int64)
            np.savetxt(handle, block, fmt="%d")
        elif n_cells:
            block = np.full((n_cells, 1 + max_n), -1, np.int64)
            block[:, 0] = sizes
            for i, element in enumerate(mesh.elements):
                block[i, 1:1 + len(element)] = element
            lines = [
                " ".join(str(int(v)) for v in row[:1 + int(row[0])])
                for row in block
            ]
            handle.write("\n".join(lines) + "\n")
        handle.write(f"CELL_TYPES {n_cells}\n")
        cell_types = np.where(sizes == 3, _VTK_TRIANGLE, np.where(sizes == 4, _VTK_QUAD, _VTK_POLYGON))
        np.savetxt(handle, cell_types, fmt="%d")
        handle.write(f"CELL_DATA {n_cells}\n")
        for name, values in (("thickness", mesh.thickness), ("eccentricity", mesh.eccentricity)):
            _scalars(handle, name, "double", values)
        _scalars(handle, "surface", "int", _codes(SURFACES, mesh.surface))
        _scalars(handle, "partition", "int", _codes(partitions, mesh.partition))
        if len(mesh.constructor) == n_cells:
            _scalars(handle, "constructor", "int", _codes(CONSTRUCTORS, mesh.constructor))
            _scalars(handle, "bay", "int", mesh.bay)
            _scalars(handle, "z0", "double", mesh.z0)
            _scalars(handle, "z1", "double", mesh.z1)
            _scalars(handle, "n0", "int", mesh.n0)
            _scalars(handle, "n1", "int", mesh.n1)
            _scalars(handle, "region", "int", _codes(regions, mesh.name))
            _scalars(handle, "flipped", "int", mesh.flipped)
        quality = getattr(mesh, "element_quality", None)
        if quality is not None and len(quality) == n_cells:
            for name in quality.dtype.names:
                _scalars(handle, name, "double", quality[name])
        hits = run_mesh_health(mesh)
        health_code = element_health_codes(hits, n_cells)
        _scalars(handle, "health_code", "int", health_code)
        sev_lookup = {name: i for i, name in enumerate(HEALTH_SEVERITIES)}
        health_sev = np.full(n_cells, len(HEALTH_SEVERITIES), np.int64)
        for hit in hits:
            if hit.element < 0 or hit.element >= n_cells:
                continue
            if health_sev[hit.element] == len(HEALTH_SEVERITIES):
                health_sev[hit.element] = sev_lookup.get(hit.severity, len(HEALTH_SEVERITIES))
        _scalars(handle, "health_severity", "int", health_sev)
        mesh.__dict__["_health_hits"] = hits
        if len(mesh.node_origin) == n_points:
            handle.write(f"POINT_DATA {n_points}\n")
            _scalars(handle, "origin", "int", _codes(ORIGINS, mesh.node_origin))
            dist = np.asarray(mesh.section_distance, float)
            if len(dist) != n_points:
                dist = np.zeros(n_points)
            _scalars(handle, "section_distance", "double", dist)
            source = np.zeros(n_points, np.int64)
            source[:int(mesh.n_source)] = 1
            _scalars(handle, "source", "int", source)
            valence = mesh.valence if len(mesh.valence) == n_points else np.zeros(n_points, np.int64)
            _scalars(handle, "valence", "int", valence)
    path.with_name(path.stem + ".legend.txt").write_text(legend + "\n", encoding="utf-8")
    if len(mesh.edges):
        _write_edges(mesh, path.with_name(path.stem + ".edges.vtk"), head)
    if len(mesh.corners):
        _write_corners(mesh, path.with_name(path.stem + ".corners.vtk"), head)
    return partitions


def _write_edges(mesh: ShellMesh, path: Path, title: str) -> None:
    edges = np.asarray(mesh.edges)
    n_edges = len(edges)
    n_points = len(mesh.nodes)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write((title + " edges")[:255] + "\n")
        handle.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        handle.write(f"POINTS {n_points} double\n")
        np.savetxt(handle, np.asarray(mesh.nodes, float), fmt="%.16g")
        handle.write(f"CELLS {n_edges} {n_edges * 3}\n")
        np.savetxt(handle, np.column_stack([np.full(n_edges, 2, np.int64), edges["a"], edges["b"]]), fmt="%d")
        handle.write(f"CELL_TYPES {n_edges}\n")
        np.savetxt(handle, np.full(n_edges, _VTK_LINE, np.int64), fmt="%d")
        handle.write(f"CELL_DATA {n_edges}\n")
        _scalars(handle, "use", "int", edges["use"])
        _scalars(handle, "length", "double", edges["length"])
        _scalars(handle, "dx", "double", edges["dx"])
        _scalars(handle, "dy", "double", edges["dy"])
        _scalars(handle, "dz", "double", edges["dz"])
        _scalars(handle, "dihedral_deg", "double", edges["dihedral_deg"])
        _scalars(handle, "sagitta", "double", edges["sagitta"])
        _scalars(handle, "kind", "int", _codes(KINDS, edges["kind"]))
        _scalars(handle, "meet", "int", _codes(MEETS, edges["meet"]))
        hits = getattr(mesh, "_health_hits", None) or run_mesh_health(mesh)
        _scalars(handle, "health_code", "int", edge_health_codes(hits, n_edges))


def _write_corners(mesh: ShellMesh, path: Path, title: str) -> None:
    """One ``VTK_VERTEX`` per corner, on its own point, so a shared node stays two vertices."""
    corners = np.asarray(mesh.corners)
    n = len(corners)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write((title + " corners")[:255] + "\n")
        handle.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        handle.write(f"POINTS {n} double\n")
        np.savetxt(handle, np.column_stack([corners["x"], corners["y"], corners["z"]]), fmt="%.16g")
        handle.write(f"CELLS {n} {n * 2}\n")
        np.savetxt(handle, np.column_stack([np.ones(n, np.int64), np.arange(n)]), fmt="%d")
        handle.write(f"CELL_TYPES {n}\n")
        np.savetxt(handle, np.full(n, _VTK_VERTEX, np.int64), fmt="%d")
        handle.write(f"CELL_DATA {n}\n")
        _scalars(handle, "element_id", "int", corners["element_id"])
        _scalars(handle, "corner", "int", corners["corner"])
        _scalars(handle, "node_id", "int", corners["node_id"])
        _scalars(handle, "angle_deg", "double", corners["angle_deg"])
        _scalars(handle, "yaw_deg", "double", corners["yaw_deg"])
        _scalars(handle, "pitch_deg", "double", corners["pitch_deg"])
        _scalars(handle, "roll_deg", "double", corners["roll_deg"])
        _scalars(handle, "bay", "int", corners["bay"])
        _scalars(handle, "constructor", "int", _codes(CONSTRUCTORS, corners["constructor"]))


__all__ = ["write_vtk"]
