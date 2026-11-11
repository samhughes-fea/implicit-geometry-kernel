"""OML and mid-surface triangle/OBJ writers for :class:`BladeGeometry`."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from implicit_geometry_kernel.blade.stations import TransverseMember
    from implicit_geometry_kernel.blade.topology.template import ShellTopologyTemplate


def write_stl_solid(fh, name: str, verts: np.ndarray, faces: list) -> None:
    fh.write(f"solid {name}\n")
    for i, j, k in faces:
        n = np.cross(verts[j] - verts[i], verts[k] - verts[i])
        ln = np.linalg.norm(n)
        if ln:
            n = n / ln
        fh.write(f" facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n  outer loop\n")
        for p in (verts[i], verts[j], verts[k]):
            fh.write(f"   vertex {p[0]:.6e} {p[1]:.6e} {p[2]:.6e}\n")
        fh.write("  endloop\n endfacet\n")
    fh.write(f"endsolid {name}\n")


def write_tri_mesh(path, verts: np.ndarray, faces: list, solid_name: str = "oml") -> None:
    path = Path(path)
    if path.suffix.lower() == ".obj":
        lines = [f"v {x:g} {y:g} {z:g}" for x, y, z in verts]
        lines += [f"f {i + 1} {j + 1} {k + 1}" for i, j, k in faces]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8") as fh:
        write_stl_solid(fh, solid_name, verts, faces)


def write_oml(geo, path, n_span: int = 41, n_chord: int = 81) -> str:
    """Watertight OML as STL or OBJ from a quad grid of span rings."""
    path = Path(path)
    zs = np.linspace(geo.stations[0].z, geo.stations[-1].z, int(n_span))
    rings = []
    for z in zs:
        poly = geo.oml_polyline(z)
        if np.allclose(poly[0], poly[-1]):
            poly = poly[:-1]
        idx = np.linspace(0, len(poly) - 1, int(n_chord), endpoint=False)
        rings.append(np.vstack([poly[int(i) % len(poly)] for i in np.round(idx)]))
    verts = np.vstack(rings)
    ns, nc = len(rings), int(n_chord)
    faces = []
    for i in range(ns - 1):
        for j in range(nc):
            a = i * nc + j
            b = i * nc + (j + 1) % nc
            c = (i + 1) * nc + (j + 1) % nc
            d = (i + 1) * nc + j
            faces.append((a, b, c))
            faces.append((a, c, d))
    for base, sign in ((0, 1), ((ns - 1) * nc, -1)):
        hub = verts[base:base + nc].mean(axis=0)
        verts = np.vstack([verts, hub])
        h = len(verts) - 1
        for j in range(nc):
            if sign > 0:
                faces.append((h, base + j, base + (j + 1) % nc))
            else:
                faces.append((h, base + (j + 1) % nc, base + j))
    write_tri_mesh(path, verts, faces, solid_name="oml")
    return str(path)


def write_mid_surfaces(
    geo,
    path,
    n_span: int = 41,
    n_chord: int = 81,
    members=None,
    include_solids: bool = False,
) -> str:
    """One lofted patch per authored region / web, as STL or OBJ (groups)."""
    from implicit_geometry_kernel.blade.topology import (
        _template_from_cut,
        loft_mapped_strips,
        rib_cell_polygons_from_samples,
        sample_region_polylines,
        section_cut_near,
    )

    path = Path(path)
    if path.suffix.lower() in {".step", ".stp"}:
        from implicit_geometry_kernel.blade.cad import write_mid_surface_step
        from implicit_geometry_kernel.blade.loft.locked_loft import locked_mid_surface_patches
        from implicit_geometry_kernel.contract import InfeasibleSectionGeometry

        z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
        zs: list[float] = [float(z) for z in np.linspace(z0, z1, int(n_span))]
        for z in geo.forced_stations():
            if z0 - 1e-12 <= float(z) <= z1 + 1e-12:
                zs.append(float(z))
        zs = sorted(set(round(float(z), 9) for z in zs))
        graphs = []
        for z in zs:
            try:
                graphs.append(geo.shell_graph(float(z), include_solids=include_solids))
            except InfeasibleSectionGeometry:
                continue
        patches = locked_mid_surface_patches(graphs, n_chord=int(n_chord), members=members)
        return write_mid_surface_step(path, [(name, rings) for name, rings, _surface in patches])

    want = None if members is None else set(members)
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = [float(z) for z in np.linspace(z0, z1, int(n_span))]
    for z in geo.forced_stations():
        if z0 - 1e-12 <= float(z) <= z1 + 1e-12:
            zs.append(float(z))
    zs = sorted(set(round(float(z), 9) for z in zs))
    samples: list[tuple[float, ShellTopologyTemplate, dict]] = []
    for z in zs:
        cut = section_cut_near(geo, float(z))
        if cut is None:
            continue
        sg = cut.shell_graph(include_solids=include_solids)
        tmpl = _template_from_cut(geo.region_graph, cut, float(z), int(n_chord))
        samples.append((float(z), tmpl, sample_region_polylines(cut, sg, tmpl, geo.region_graph)))
    groups: dict[str, list[tuple[np.ndarray, bool]]] = {}
    for (_za, t0, s0), (_zb, t1, s1) in zip(samples[:-1], samples[1:]):
        for name, face, collapsing in loft_mapped_strips(s0, s1, t0, t1):
            if want is not None and name not in want:
                continue
            groups.setdefault(name, []).append((face, collapsing))
    chunks: list[tuple[str, np.ndarray, list]] = []
    for name, faces in groups.items():
        verts: list[np.ndarray] = []
        idx: list[tuple[int, ...]] = []
        for face, _collapsing in faces:
            base = len(verts)
            for p in face:
                verts.append(np.asarray(p, float).reshape(3))
            n = len(face)
            if n == 3:
                idx.append((base, base + 1, base + 2))
            elif n >= 4:
                idx.append(tuple(base + i for i in range(n)))
        if verts:
            chunks.append((name, np.vstack(verts), idx))
    # rib faces (a "rib" member filter keeps them all) and the tip closure, which is
    # part of the skin (a "skin" or "tip_closure" filter keeps it)
    tip = getattr(geo, "tip_closure", None)
    for z, tmpl, samp in samples:
        closures: list[tuple[TransverseMember | str, tuple[str, ...]]] = [
            (m, ("rib",)) for m in geo.transverse if abs(float(m.z) - float(z)) <= 1e-8
        ]
        if tip is not None and abs(float(geo.tip_closure_z) - float(z)) <= 1e-8:
            closures.append((str(tip.name), ("tip_closure", "skin")))
        for member, filters in closures:
            for name, ring in rib_cell_polygons_from_samples(tmpl, samp, member):
                if want is not None and name not in want and not want.intersection(filters):
                    continue
                n = len(ring)
                if n < 3:
                    continue
                chunks.append((name, np.asarray(ring, float), [tuple(range(n))]))
    if path.suffix.lower() == ".obj":
        lines: list[str] = []
        base = 0
        for name, chunk_verts, chunk_faces in chunks:
            lines.append(f"g {name}")
            lines += [f"v {x:g} {y:g} {z:g}" for x, y, z in chunk_verts]
            for face in chunk_faces:
                lines.append("f " + " ".join(str(base + i + 1) for i in face))
            base += len(chunk_verts)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    else:
        with path.open("w", encoding="utf-8") as fh:
            for name, chunk_verts, chunk_faces in chunks:
                tris: list[tuple[int, ...]] = []
                for face in chunk_faces:
                    if len(face) == 3:
                        tris.append(face)
                    else:
                        tris.extend((face[0], face[i], face[i + 1]) for i in range(1, len(face) - 1))
                write_stl_solid(fh, name, chunk_verts, tris)
    return str(path)
