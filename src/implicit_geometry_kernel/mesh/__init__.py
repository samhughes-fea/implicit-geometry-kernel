"""Meshing backends. Geometry stays in the kernel.

Layout:

* ``mesh.shell`` — conforming shell mesh
* ``mesh.health`` — classify/flag registry (engineer judgement)
* ``mesh.quality`` — per-element metrics and check rollups
* ``mesh.tables`` — edges/corners/CSV dump
* ``mesh.export`` — VTK writers
* ``mesh.fill`` — polygon Voronoi and Steiner Delaunay
* ``mesh.gmsh`` — optional Gmsh comparison

``shell_mesh`` is the blade's conforming shell mesh (quads on the skin and webs,
triangles on ribs and the tip). Import from this package or the subpackages
above (``mesh.shell``, ``mesh.export``, ``mesh.fill``). Flat module paths
(``mesh.shell_mesh``, ``mesh.vtk``, ``mesh.polygon``, …) were removed in 0.9.x (K5).
"""

from implicit_geometry_kernel.mesh.export import write_vtk
from implicit_geometry_kernel.mesh.fill import polygon_mesh_ring
from implicit_geometry_kernel.mesh.health import (
    HEALTH_CLASSES, HEALTH_CODES, HealthClass, HealthHit, HealthThresholds,
    diff_study_folders, run_mesh_health,
)
from implicit_geometry_kernel.mesh.quality import MeshCheck, MeshRunQuality, mesh_check, mesh_quality
from implicit_geometry_kernel.mesh.shell import MeshRecipe, ShellMesh, mesh_stations, shell_mesh
from implicit_geometry_kernel.mesh.tables import (
    mesh_failures, off_root_open_edges, trailing_face_failures, transition_band_failures,
    write_mesh_tables,
)

__all__ = [
    "HEALTH_CLASSES", "HEALTH_CODES", "HealthClass", "HealthHit", "HealthThresholds",
    "MeshCheck", "MeshRecipe", "MeshRunQuality", "ShellMesh",
    "diff_study_folders", "mesh_check", "mesh_failures", "mesh_quality", "mesh_stations",
    "off_root_open_edges", "polygon_mesh_ring", "run_mesh_health", "shell_mesh",
    "trailing_face_failures", "transition_band_failures", "write_mesh_tables", "write_vtk",
]
