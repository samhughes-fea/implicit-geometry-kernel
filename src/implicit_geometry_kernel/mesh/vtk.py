"""Write a :class:`~implicit_geometry_kernel.mesh.shell_mesh.ShellMesh` as legacy ASCII VTK.

VTK keeps polygons as polygons (``VTK_POLYGON``), so rib and tip cells are
written as they are, with no triangle split. Cell data: ``thickness``,
``eccentricity``, and integer ids for ``surface`` and ``partition`` (their
names are listed in the file header comment order, see ``partition_names``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from implicit_geometry_kernel.mesh.shell_mesh import SURFACES, ShellMesh

_VTK_QUAD, _VTK_POLYGON = 9, 7


def write_vtk(mesh: ShellMesh, path: str | Path, title: str = "shell mesh") -> tuple[str, ...]:
    """Write ``mesh`` to ``path``; returns the partition names in id order."""
    partitions = tuple(dict.fromkeys(mesh.partition))
    part_id = {p: i for i, p in enumerate(partitions)}
    n_cells = len(mesh.elements)
    size = sum(len(e) + 1 for e in mesh.elements)
    lines = ["# vtk DataFile Version 3.0", " ".join(str(title).split())[:255] or "shell mesh", "ASCII",
             "DATASET UNSTRUCTURED_GRID", f"POINTS {len(mesh.nodes)} double"]
    lines += [f"{x!r} {y!r} {z!r}" for x, y, z in np.asarray(mesh.nodes, float).tolist()]
    lines.append(f"CELLS {n_cells} {size}")
    lines += [" ".join(str(k) for k in (len(e), *e)) for e in mesh.elements]
    lines.append(f"CELL_TYPES {n_cells}")
    lines += [str(_VTK_QUAD if len(e) == 4 else _VTK_POLYGON) for e in mesh.elements]
    lines.append(f"CELL_DATA {n_cells}")
    for name, values in (("thickness", mesh.thickness), ("eccentricity", mesh.eccentricity)):
        lines += [f"SCALARS {name} double 1", "LOOKUP_TABLE default"]
        lines += [repr(float(v)) for v in values]
    lines += ["SCALARS surface int 1", "LOOKUP_TABLE default"]
    lines += [str(SURFACES.index(s)) for s in mesh.surface]
    lines += ["SCALARS partition int 1", "LOOKUP_TABLE default"]
    lines += [str(part_id[p]) for p in mesh.partition]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return partitions


__all__ = ["write_vtk"]
