"""Meshing backends. Geometry stays in the kernel.

``polygon`` is a kernel Voronoi fill (GEOS). ``delaunay`` is a numpy Steiner
trimesh. ``gmsh`` is the optional extra.
"""

from implicit_geometry_kernel.mesh.gmsh import mesh_rib_cells, mesh_shell, write_shell_msh
from implicit_geometry_kernel.mesh.polygon import polygon_mesh_ring

__all__ = ["mesh_rib_cells", "mesh_shell", "polygon_mesh_ring", "write_shell_msh"]
