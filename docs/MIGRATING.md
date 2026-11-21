# Migrating from 0.2.x to 0.9

Code written against 0.2.x stops at the first of these changes. Each row gives
the 0.2.x call and its 0.9 form. Geometry is unchanged by the migration: on a
NACA 0024 two-web blade, web lengths, cap length and cell areas agree with
0.2.0 to 2e-4 relative (`tests/test_compat_0_2.py`).

## Install

| 0.2.x | 0.9 |
|---|---|
| Shapely 2 pulled in by pip | Install GEOS first: `apt install libgeos-dev`, `brew install geos` or `conda install geos`. The kernel loads `libgeos_c` through ctypes; set `GEOS_LIBRARY_PATH` if it is not on the loader path |

## Layout codes

A layout code names the section family, the web anchoring and the rib count.
The short form is rejected with `KeyError`.

| 0.2.x | 0.9 |
|---|---|
| `"2D-F"` | `"2D-F-LC-0R"` |
| `"0A"`, `"3B-CN"`, ... | `"0A-LC-0R"`, `"3B-CN-LC-0R"`, ... |

Anchoring is `LC` (local chord fraction), `LE` (distance from the leading
edge) or `PA` (blade-x from the pitch axis). Web positions given as chord
fractions match 0.2.x with `LC`. `0R` means no ribs; `10R` adds ten interior
ribs. Every blade also has a tip closure, which is not a rib (0.9.2).

## BladeGeometry

| 0.2.x | 0.9 |
|---|---|
| `BladeGeometry(stations, "2D-F", t_skin, t_cap, t_web, ...)` | `BladeGeometry(stations, "2D-F-LC-0R", t_skin, t_cap, t_web, ..., pivot_fraction=0.333)` |
| no `pivot_fraction` (fell back to 0.333) | `pivot_fraction` is required: a chord fraction or a span table. Use the stations' pitch axis where you set one |
| `geo.station_at(z)` | `geo.interpolate(z)`: `(chord, twist, pitch_axis, offset, (x, upper, lower))` |

## Imports

| 0.2.x | 0.9 |
|---|---|
| `implicit_geometry_kernel.section_geometry.layout.system_layout.SYSTEM_TYPE_KEYS` | `implicit_geometry_kernel.SYSTEM_TYPE_KEYS` |
| `implicit_geometry_kernel.blade_geometry` | `implicit_geometry_kernel.blade` |

## Section level

| 0.2.x | 0.9 |
|---|---|
| `SectionGeometry(layout=..., airfoil=..., ...)` | `SectionGeometry(SectionSpec.from_kwargs(...), SectionOptions.from_kwargs(...))` |
| `apply_shapely` | `apply_clip` |
| `cut.mid_surface()` put a cap + skin run on the stack's mid-surface | The run stays on the skin mid-surface. Give the stacked laminate an offset of half the cap thickness, inward (`cut.shell_graph().offset_regions()`) |

Changes made after 0.2 and removed again before 0.9 (`span_kind`, `SpanWeb`,
`ExprSDF`, `WallRegionGraph.cap_mode`) are listed in the changelog.

## 0.9.1 to 0.9.2: the tip closure

Every `BladeGeometry` now closes its skin at the tip with a `TipClosure`. It is
part of the skin, not a rib.

| 0.9.1 | 0.9.2 |
|---|---|
| `-10R` gave 9 ribs plus a `tip_closure` rib at the tip | 10 interior ribs at `z0 + k (z1 - z0) / 11`, plus the tip closure |
| `TransverseMember("tip_closure", z_tip, t)` | `BladeGeometry(..., tip_closure=TipClosure(thickness=t))`. The old member still works at the tip, with a `DeprecationWarning` |
| tip closure thickness 0.012 m (catalogue) or the windIO skin/TE gauge | the skin thickness law at the tip, unless `TipClosure(thickness=...)` sets it |
| `transverse_members_from_windio` returned a `tip_closure` member first | it returns the YAML ribs/bulkheads only |
| tip face in `shell_graph().faces` with `surface="rib"` | `surface="skin"`, `kind="tip_closure"`, `layers=("skin",)` |
| `rib_face_polygons(geo)` included the tip | ribs only; `tip_closure_polygons(geo)` gives the tip closure |

## 0.9.1 to 0.9.2: the shell mesh

| 0.9.1 | 0.9.2 |
|---|---|
| `mesh.gmsh.mesh_shell` / `write_shell_msh` (quads plus triangles, own nodes per strip) | `mesh.shell_mesh(geo, n_span \| h=)`: one conforming sheet, all-quad skin and webs, polygon ribs and tip; `mesh.write_vtk` to write it. The Gmsh calls warn `DeprecationWarning` |
| `mesh.gmsh.triangulate_ring`, `fill_rib_cells`, `mesh_rib_cells` (triangles) | `mesh.polygon_mesh_ring` (cells of 5 or more nodes), or the rib elements of `shell_mesh` |
| `sample_region_polylines` gave `n_chord` nodes per region | the next odd count: `n_chord + 1` when `n_chord` is even (TE regions and the wrap round up likewise) |
| mesh `cap_*` regions bounded by one chord window on both sides | bounded by each side's cap edges in the cut; on a twisted blade the region keys move by up to ~0.1 chord |
| a cap whose skin is concave under it could be missing from `shell_graph().loop` | it is always there (SG-029) |

## 0.9.2 to 0.9.3: flat mesh modules removed (K5)

Import the conforming shell mesh and writers from the package or subpackages.
0.9.3 removed the old flat module paths; 0.9.4 brings them back with a
`DeprecationWarning` on each name read, and 0.9.5 removes them.

| Removed module path | Use instead |
|---|---|
| `implicit_geometry_kernel.mesh.shell_mesh` | `from implicit_geometry_kernel.mesh import shell_mesh` or `mesh.shell` |
| `implicit_geometry_kernel.mesh.vtk` | `from implicit_geometry_kernel.mesh import write_vtk` or `mesh.export` |
| `implicit_geometry_kernel.mesh.polygon` / `polygon_voronoi` / `delaunay` | `mesh.fill` (`polygon_mesh_ring`, `triangulate_ring`, …) |
| `mesh.mesh_rib_cells`, `mesh.mesh_shell`, `mesh.write_shell_msh` at the package root | `mesh.gmsh` (comparison only) |
| `ShellMesh` per-element fields as lists | tuples; convert with `list(...)` where a list is needed |
| rib and tip-closure cells: one polygon of 5 or more nodes each | triangles (`len(element) == 3`); a consumer that fanned polygons takes them as they are |

Maintainer reports live under `docs/internal/` (stubs remain at the old paths).

## 0.9.2 to 0.9.3: top-level barrel trim (K5)

Rare topology helpers left the top-level `__all__`. Import them from
`implicit_geometry_kernel.blade` (or `blade.topology`). A root import still
works in 0.9.x with a `DeprecationWarning`.

| Left the barrel | Import from |
|---|---|
| `TopologyDiff`, `TopologyEvent`, `topology_diff` | `blade` |
| `ShellTopologyTemplate`, `Region` | `blade` / `blade.topology` |
| `PartitionEdge`, `partition_edge_curve`, `ShellRun`, `MidSurface` | `blade` |
| `CellPolygon`, `SectionWall` | `blade` |
| `shared_loft_names`, `signatures_compatible`, `zone_region_graphs` | `blade` |
| `shell_strips_from_graph`, `station_network_from_meta`, `refine_section_samples` | `blade` |

Loft, station, SDF and catalogue names the consumer already imports at the
root stay on `__all__`.

## Checking a consumer

Run the consumer's tests against a kernel checkout:

```bash
python scripts/check_consumer.py PATH/TO/CONSUMER
```
