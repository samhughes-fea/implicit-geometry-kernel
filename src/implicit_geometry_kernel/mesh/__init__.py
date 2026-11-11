"""Meshing backends. Geometry stays in the kernel.

``shell_mesh`` is the blade's conforming shell mesh (quads on the skin and webs,
polygons on ribs and the tip). ``polygon`` is a kernel Voronoi fill (GEOS). ``delaunay`` is a numpy Steiner
trimesh. ``gmsh`` is the optional extra.
"""

from implicit_geometry_kernel.mesh.gmsh import mesh_rib_cells, mesh_shell, write_shell_msh
from implicit_geometry_kernel.mesh.polygon import polygon_mesh_ring
from implicit_geometry_kernel.mesh.shell_mesh import ShellMesh, mesh_stations, shell_mesh

__all__ = ["ShellMesh", "mesh_rib_cells", "mesh_shell", "mesh_stations", "polygon_mesh_ring", "shell_mesh",
           "write_shell_msh"]
