"""Gmsh backend (optional extra ``[mesh]``)."""

from implicit_geometry_kernel.mesh.gmsh.plane import (
    fill_rib_cells,
    mesh_rib_cells,
    mesh_rib_cells_many,
    triangulate_ring,
)
from implicit_geometry_kernel.mesh.gmsh.shell import (
    junction_web_errors,
    mesh_shell,
    write_shell_msh,
)

__all__ = [
    "fill_rib_cells",
    "junction_web_errors",
    "mesh_rib_cells",
    "mesh_rib_cells_many",
    "mesh_shell",
    "triangulate_ring",
    "write_shell_msh",
]
