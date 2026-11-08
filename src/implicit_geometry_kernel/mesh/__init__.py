"""Optional shell meshing (Gmsh). Geometry stays in the kernel; this is the backend."""

from implicit_geometry_kernel.mesh.gmsh_shell import mesh_shell, write_shell_msh

__all__ = ["mesh_shell", "write_shell_msh"]
