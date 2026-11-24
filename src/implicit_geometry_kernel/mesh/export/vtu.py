"""Write a :class:`~implicit_geometry_kernel.mesh.shell.ShellMesh` as binary XML VTK (``.vtu``).

The file is a VTK ``UnstructuredGrid`` with base64-encoded binary data arrays
(``format="binary"``, a ``UInt64`` byte-count header per array, little-endian),
which ParaView and the ``vtk`` package read directly. It is several times
smaller and faster to write than the legacy ASCII file from
:func:`~implicit_geometry_kernel.mesh.export.write_vtk`.

Cell data: ``thickness``, ``eccentricity`` and integer codes for ``surface``,
``partition`` and ``constructor``. The code legends and the mesh
:meth:`~implicit_geometry_kernel.mesh.shell.ShellMesh.provenance` are written as
JSON in an XML comment at the top of the file, so the file says which mesh it is.
"""
from __future__ import annotations

import base64
import json
from pathlib import Path

import numpy as np

from implicit_geometry_kernel.mesh.shell import SURFACES, ShellMesh
from implicit_geometry_kernel.mesh.tables import CONSTRUCTORS

_VTK_TRIANGLE, _VTK_QUAD, _VTK_POLYGON = 5, 9, 7
_TYPES = {np.dtype("<f8"): "Float64", np.dtype("<i8"): "Int64", np.dtype("u1"): "UInt8"}


def _codes(names, values) -> np.ndarray:
    lookup = {name: i for i, name in enumerate(names)}
    return np.fromiter((lookup.get(str(v), len(names)) for v in values), np.int64, count=len(values))


def _array(name: str, values: np.ndarray, components: int = 1) -> str:
    arr = np.ascontiguousarray(values)
    raw = arr.tobytes()
    header = np.asarray([len(raw)], dtype="<u8").tobytes()
    payload = base64.b64encode(header + raw).decode("ascii")
    kind = _TYPES[arr.dtype]
    comp = f' NumberOfComponents="{components}"' if components != 1 else ""
    return f'<DataArray type="{kind}" Name="{name}"{comp} format="binary">{payload}</DataArray>'


def write_vtu(mesh: ShellMesh, path: str | Path) -> tuple[str, ...]:
    """Write ``mesh`` to ``path`` as binary ``.vtu``; returns the partition names in id order."""
    path = Path(path)
    partitions = tuple(dict.fromkeys(mesh.partition))
    n_points, n_cells = len(mesh.nodes), len(mesh.elements)
    sizes = np.fromiter((len(e) for e in mesh.elements), np.int64, count=n_cells)
    connectivity = np.fromiter((i for e in mesh.elements for i in e), np.int64, count=int(sizes.sum()))
    offsets = np.cumsum(sizes).astype("<i8")
    types = np.where(sizes == 3, _VTK_TRIANGLE, np.where(sizes == 4, _VTK_QUAD, _VTK_POLYGON)).astype("u1")
    meta = {
        "provenance": mesh.provenance(),
        "surface": list(SURFACES),
        "partition": list(partitions),
        "constructor": list(CONSTRUCTORS),
    }
    comment = json.dumps(meta, separators=(",", ":")).replace("--", "- -")
    cell_data = [
        _array("thickness", np.asarray(mesh.thickness, "<f8")),
        _array("eccentricity", np.asarray(mesh.eccentricity, "<f8")),
        _array("surface", _codes(SURFACES, mesh.surface).astype("<i8")),
        _array("partition", _codes(partitions, mesh.partition).astype("<i8")),
    ]
    if len(mesh.constructor) == n_cells:
        cell_data.append(_array("constructor", _codes(CONSTRUCTORS, mesh.constructor).astype("<i8")))
    parts = [
        '<?xml version="1.0"?>',
        f"<!-- {comment} -->",
        '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">',
        "<UnstructuredGrid>",
        f'<Piece NumberOfPoints="{n_points}" NumberOfCells="{n_cells}">',
        "<Points>",
        _array("Points", np.asarray(mesh.nodes, "<f8").reshape(-1, 3), components=3),
        "</Points>",
        "<Cells>",
        _array("connectivity", connectivity.astype("<i8")),
        _array("offsets", offsets),
        _array("types", types),
        "</Cells>",
        "<CellData>",
        *cell_data,
        "</CellData>",
        "</Piece>",
        "</UnstructuredGrid>",
        "</VTKFile>",
    ]
    path.write_text("\n".join(parts) + "\n", encoding="ascii")
    return partitions


def read_vtu_arrays(path: str | Path) -> dict:
    """Read back a file from :func:`write_vtu`: ``{name: array}`` plus ``"_meta"`` (the JSON comment).

    A small reader for tests and round trips; it handles only this writer's layout.
    """
    import re

    text = Path(path).read_text(encoding="ascii")
    out: dict = {}
    m = re.search(r"<!-- (.*?) -->", text, re.S)
    if m:
        out["_meta"] = json.loads(m.group(1).replace("- -", "--"))
    kinds = {"Float64": "<f8", "Int64": "<i8", "UInt8": "u1"}
    for kind, name, comp, payload in re.findall(
        r'<DataArray type="(\w+)" Name="(\w+)"(?: NumberOfComponents="(\d+)")? format="binary">([^<]*)</DataArray>',
        text,
    ):
        raw = base64.b64decode(payload)
        n = int(np.frombuffer(raw[:8], "<u8")[0])
        arr = np.frombuffer(raw[8:8 + n], kinds[kind])
        out[name] = arr.reshape(-1, int(comp)) if comp else arr
    return out


__all__ = ["read_vtu_arrays", "write_vtu"]
