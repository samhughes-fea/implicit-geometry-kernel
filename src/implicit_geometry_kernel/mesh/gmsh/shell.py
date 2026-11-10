"""Gmsh wall-strip mesh from authored templates and feature-curve vertices."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from implicit_geometry_kernel.blade.topology import (
    loft_mapped_strips,
    sample_region_polylines,
    section_cut_near,
    templates_along,
)
from implicit_geometry_kernel.extract.curves import locate_shell_vertex


def _stations(geo, zs, n_chord: int):
    tpls = templates_along(geo, zs, n_chord=n_chord)
    walls = []
    for tmpl in tpls:
        cut = section_cut_near(geo, tmpl.z)
        if cut is None:
            continue
        samples = sample_region_polylines(cut, cut.shell_graph(), tmpl, geo.region_graph)
        if samples:
            walls.append((float(tmpl.z), tmpl, samples, cut))
    return walls


def _nid(p, index: dict, coords: list) -> int:
    key = tuple(round(float(c), 8) for c in np.asarray(p, float).ravel()[:3])
    if key not in index:
        index[key] = len(index) + 1
        coords.append(key)
    return index[key]


def _unique_ring(tags: list[int]) -> list[int]:
    out: list[int] = []
    for t in tags:
        if not out or t != out[-1]:
            out.append(t)
    if len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out


def mesh_shell(geo, zs=None, *, n_chord: int = 8, out: Path | None = None) -> dict:
    """Build a Gmsh shell of mapped strips. Returns node/element counts.

    Compatible bands become quads (one per authored strip). Birth/death faces
    stay triangles until all-quad subdivision. Requires the ``gmsh`` package.
    """
    import gmsh

    if zs is None:
        z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
        zs = np.linspace(z0, z1, 21)
    walls = _stations(geo, zs, n_chord)
    if len(walls) < 2:
        raise ValueError("need at least two feasible stations to mesh")
    index: dict[tuple, int] = {}
    coords: list[tuple] = []
    tris: list[list[int]] = []
    quads: list[list[int]] = []
    for (_z0, t0, s0, _c0), (_z1, t1, s1, _c1) in zip(walls, walls[1:]):
        for _name, face, collapsing in loft_mapped_strips(s0, s1, t0, t1):
            tags = _unique_ring([_nid(p, index, coords) for p in np.asarray(face, float)])
            if collapsing or len(tags) == 3:
                if len(tags) >= 3:
                    tris.append(tags[:3])
            elif len(tags) == 4:
                quads.append(tags)
            elif len(tags) > 4:
                for i in range(1, len(tags) - 1):
                    tris.append([tags[0], tags[i], tags[i + 1]])
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("shell")
        gmsh.model.addDiscreteEntity(2, 1)
        tags = list(range(1, len(coords) + 1))
        flat = [c for xyz in coords for c in xyz]
        gmsh.model.mesh.addNodes(2, 1, tags, flat)
        if tris:
            gmsh.model.mesh.addElementsByType(1, 2, [], [n for t in tris for n in t])
        if quads:
            gmsh.model.mesh.addElementsByType(1, 3, [], [n for q in quads for n in q])
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
        if tris or quads:
            gmsh.model.mesh.refine()
        n_nodes = gmsh.model.mesh.getNodes()[0].size
        n_elem = int(sum(len(t) for t in gmsh.model.mesh.getElements(2)[1]))
        faces = _element_faces(gmsh)
        if out is not None:
            Path(out).parent.mkdir(parents=True, exist_ok=True)
            gmsh.write(str(out))
        return {
            "n_nodes": int(n_nodes),
            "n_elements": n_elem,
            "n_stations": len(walls),
            "faces": faces,
        }
    finally:
        gmsh.finalize()


def _element_faces(gmsh) -> list[np.ndarray]:
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = np.asarray(coords, float).reshape(-1, 3)
    lookup = {int(t): xyz[i] for i, t in enumerate(tags)}
    n_per = {2: 3, 3: 4}
    faces = []
    types, _, node_tags = gmsh.model.mesh.getElements(2)
    for etype, nodes in zip(types, node_tags):
        k = n_per.get(int(etype))
        if k is None:
            continue
        for row in np.asarray(nodes, int).reshape(-1, k):
            faces.append(np.array([lookup[int(n)] for n in row]))
    return faces


def write_shell_msh(geo, path: Path, zs=None, *, n_chord: int = 8) -> dict:
    return mesh_shell(geo, zs, n_chord=n_chord, out=Path(path))


def junction_web_errors(geo, zs, *, n_chord: int = 8) -> float:
    """Max distance from a web-foot vertex to the sampled web end, / chord."""
    walls = _stations(geo, zs, n_chord)
    worst = 0.0
    for z, tmpl, samples, cut in walls:
        chord = max(float(cut.chord), 1e-9)
        for web in tmpl.webs:
            if web.name not in samples:
                continue
            ends = np.asarray(samples[web.name], float)
            for vid, end in ((web.suction, ends[0]), (web.pressure, ends[-1])):
                loc = locate_shell_vertex(geo, vid, [z], n_chord=n_chord)
                if len(loc):
                    worst = max(worst, float(np.linalg.norm(loc[0] - end[:3])) / chord)
    return worst
