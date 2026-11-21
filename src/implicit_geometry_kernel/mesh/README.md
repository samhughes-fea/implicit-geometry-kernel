# mesh

Conforming shell meshes, fills, quality metrics, and health classification.

↑ [package](../README.md) · [docs hub](../../../docs/README.md)

## Layout

| Subpackage | Role |
|---|---|
| `shell/` | Blade shell mesh (`shell_mesh`, stations, recipe) |
| `health/` | Flag registry (`run_mesh_health`, `HealthHit`) — classify only |
| `quality/` | Element metrics and `MeshCheck` rollups |
| `tables/` | Edges, corners, seams, CSV dump |
| `export/` | VTK writers |
| `fill/` | Polygon Voronoi and Steiner Delaunay |
| `gmsh/` | Optional Gmsh comparison (`[mesh]` extra) |

## Imports

```python
from implicit_geometry_kernel.mesh import shell_mesh, write_mesh_tables, run_mesh_health
from implicit_geometry_kernel.mesh import HEALTH_CLASSES, diff_study_folders
from implicit_geometry_kernel.mesh.shell import ShellMesh
from implicit_geometry_kernel.mesh.health import HealthHit
from implicit_geometry_kernel.mesh.export import write_vtk
from implicit_geometry_kernel.mesh.fill import polygon_mesh_ring, triangulate_ring
```

Prefer the subpackages above. Flat module paths (`mesh.shell_mesh`, `mesh.vtk`,
`mesh.polygon`, `mesh.delaunay`, `mesh.polygon_voronoi`) were removed; use
`from implicit_geometry_kernel.mesh import shell_mesh, write_vtk, …` or
`mesh.shell` / `mesh.export` / `mesh.fill`. `run_mesh_health` flags only; the engineer judges the dump. See the
[glossary](../../../docs/glossary.md) for registry codes.
