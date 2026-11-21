# Glossary

↑ [docs](README.md)

Look-up of public types. One home page per name. Barrel names are those in
`implicit_geometry_kernel.__all__`. Station and contract types are imported from
their packages. Rare topology helpers (`TopologyDiff`, `ShellTopologyTemplate`,
`PartitionEdge`, …) left the barrel in 0.9.x (K5); their home remains
[blade](blade/README.md) / `blade.topology`.

Do not merge lookalike names.

## Contents

1. [Loft](#loft)
2. [Station](#station)
3. [Contract](#contract)
4. [SDF](#sdf)
5. [Mesh health](#mesh-health)
6. [Lookalikes](#lookalikes)
7. [Barrel functions](#barrel-functions)

---

## Loft

Home: [blade](blade/README.md#in-memory-objects). Source: `blade/`.

| Name | What it is | File |
|---|---|---|
| `Airfoil` | Aerofoil vertices or family (NACA, CST, ellipse) | `blade/airfoil.py` |
| `Station` | Span \(z\), chord, twist, pitch axis, offset, aerofoil | `blade/stations.py` |
| `SkinRegion` | Chordwise skin band on both surfaces | `blade/stations.py` |
| `ShearWeb` | Web at a chord fraction (station DTO) | `blade/stations.py` |
| `Stringer` | Partial-span or free-end branch | `blade/stations.py` |
| `TransverseMember` | Rib, bulkhead, or frame at one \(z\) | `blade/stations.py` |
| `PlyDrop` | Geometry-only extra thickness window in \((z,s)\) | `blade/stations.py` |
| `BladeGeometry` | Loft interpolants; not a snapshot of every wall | `blade/loft/loft.py` |
| `section_at` | Function wrapping `BladeGeometry.section_at` | `blade/loft/loft.py` |
| `SectionCut` | Frozen cut at \(z\): walls, cells, dropped webs | `blade/cut.py` |
| `station_webs` | Realised web count W on a `SectionCut` after drop | `blade/cut.py` |
| `station_family` | Realised cap family X on a `SectionCut` after drop | `blade/cut.py` |
| `SectionWall` | One wall midline in the blade frame (blade package; not on barrel) | `blade/cut.py` |
| `CellPolygon` | Closed cell interior on a loft cut (blade package; not on barrel) | `blade/cut.py` |
| `ShellGraph` | Loftable skin loop, web branches, rib faces | `blade/shell_graph.py` |
| `ShellRun` | One loftable run on the graph (blade package; not on barrel) | `blade/shell_graph.py` |
| `MidSurface` | Oriented mid-surface segments of a cut (blade package; not on barrel) | `blade/shell_graph.py` |
| `PartitionEdge` | Meeting vertex between two loop partitions (blade package; not on barrel) | `blade/shell_graph.py` |
| `shell_graph` | Builder: `SectionCut` → `ShellGraph` | `blade/shell_graph.py` |
| `BladeConcept` | Catalogue or parametric layout plus spanwise extras | `blade/concepts/concept.py` |
| `SpanwiseScheme` | Ribs, stringers, ply drops compiled to kwargs | `blade/concepts/concept.py` |
| `LockedLoftPlan` | Held point count per run across topology events | `blade/loft/locked_loft.py` |
| `resample_locked` | Resample graph runs to a locked plan | `blade/loft/locked_loft.py` |
| `shell_loft_zones` | Name-signature bands along the span | `blade/loft/locked_loft.py` |
| `StackThickness` | Nearest-run stack \(T\) on `shell_graph`, ramped at cap edges | `blade/topology/fields.py` |
| `CutWallMetaSDF` | Spanwise `wall_meta` DAG (`section_at` per unique \(z\)) | `blade/topology/fields.py` |
| `WallRegionGraph` | Authored loop regions and webs | `blade/topology/graph.py` |
| `Region` | Named loop partition on the topology graph (blade package; not on barrel) | `blade/topology/graph.py` |
| `Anchor` | Named vertex on a region graph | `blade/topology/graph.py` |
| `WebMember` | Web on the region graph | `blade/topology/graph.py` |
| `ShellTopologyTemplate` | Authored topology template (blade package; not on barrel) | `blade/topology/template.py` |
| `StationNetwork` | Merged wall graph from feature meta | `blade/topology/network.py` |
| `TopologyEvent` | Onset or drop of a named run (blade package; not on barrel) | `blade/shell_graph.py` |
| `TopologyDiff` | Difference of two graph signatures (blade package; not on barrel) | `blade/shell_graph.py` |

## Station

Home: [catalogue](catalogue/README.md) for tokens; [section](section/README.md) for the worker.
Not all are on the kernel barrel.

| Name | What it is | File |
|---|---|---|
| `SectionGeometry` | One-station worker | `section/spec/` |
| `SectionGeometryResult` | Gated products (`properties_map`, midlines, loops) | `section/spec/results.py` |
| `SystemLayoutSpec` | Frozen catalogue token | `catalogue/system_layout.py` |
| `SectionLayout` | Parametric layout (or catalogue preset) | `section/layout/section_layout.py` |
| `CapPlacementKind` | Cap plan: `none` / `fixed_band` / `per_web_pads` / `continuous_box` | `section/layout/cap_placement.py` |
| `CapPlacementRule` | Frozen cap kind plus window / pad fields | `section/layout/cap_placement.py` |
| `SYSTEM_TYPE_KEYS` | All registered WX-Y-Z-NR tokens | `catalogue/system_layout.py` |
| `SYSTEM_TYPE_LAYOUT_KEYS` | Forty `*-LC-0R` cores for section goldens | `catalogue/system_layout.py` |
| `resolve_system_type` | Token → `SystemLayoutSpec` | `catalogue/system_layout.py` |
| `MultiCellSection` | Assembled named solids (`SectionMidlineHost`) | `section/assembly/` |
| `StationLayoutResolution` | Effective layout after web drops | `section/spec/` |

## Contract

Home: [contract](contract/README.md). Not on the kernel barrel.

| Name | What it is | File |
|---|---|---|
| `ShellMidlineStrip` | Mid-surface polyline in \(\mathrm{B}\) plus thickness | `contract/shell_midline_strip.py` |
| `StripBond` | Bond between two strips | `contract/shell_midline_strip.py` |
| `CellWallLoop` | Closed cell circuit for the beam/section model | `contract/cell_wall_loop.py` |
| `LaminateOffsetRegion` | Constitutive OFFSET (cap eccentricity), not a kinked polyline | `contract/laminate_offset_region.py` |
| `InfeasibleSectionGeometry` | Unbuildable station (cavity too thin, last web dropped, …) | `contract/errors.py` |

## SDF

Home: [sdf](sdf/README.md). Barrel names:

| Name | What it is |
|---|---|
| `AirfoilSDF` | Implicit polygon field of an aerofoil contour (`sdf.shapes`) |
| `SDFGrid` | Evaluation grid (sample / gradient / debug contour) |
| `SDFNode` | Node of the CSG expression tree (`sdf.nodes`) |
| `FieldSample` | Sampled field plus optional feature meta |
| `FeatureTable` | Integer ids for tagged walls |
| `NO_FEATURE` | Sentinel id |
| `TagNode` | IR node that stamps a feature id |
| `Circle` / `Polygon` | Canonical 2-D Expr leaves (`SDFNode`) |
| `PriorityUnionNode` | Union resolved by feature priority |
| `KIND_PRIORITY` | Default meta ownership by kind (cap beats skin) on `FeatureTable` |
| `LoftedPolygonSDF` | Polygon SDF at each point's `z`; `accelerate` BVH via `accel_at` / `bvh` |
| `GeometrySolver` | In-plane Newton (`sdf.solve`); not the extract package |

## Mesh health

Home: [`mesh/README.md`](../src/implicit_geometry_kernel/mesh/README.md). Source: `mesh/health/`.
Flags **classify** observations; they do not pass or fail a mesh.

| Name | What it is | Inspect next |
|---|---|---|
| `HealthClass` / `HEALTH_CLASSES` | Registry entry: code, layer, severity hint, scope, doc | `HEALTH_CODES` order matches VTK `health_code` |
| `HealthHit` | One flag on an element, edge, or node | `failures.csv` |
| `HealthThresholds` | Cut-offs used only to classify | kwargs to `run_mesh_health` |
| `run_mesh_health` | Emit all registered hits for a shell mesh | `failures.csv`, `health_counts.csv` |
| `diff_study_folders` | Before/after count and element-id delta | two example-15 figure folders |
| `open_edge_off_root` | Free edge off the root station | `edges.csv` + region/kind/z |
| `t_junction` / `near_miss` | Node on a foreign edge / coincident pair | `junctions.csv` |
| `quad_angle` / `quad_inverted` | Shallow or inverted quad | `elements.csv` angle_min / jacobian_min |
| `trailing_face_*` | Region `te` wrap quality | `elements.csv` region=te |
| `sheet_intersection` | Non-adjacent sheets pierce | VTK both regions |
| `ridge_mismatch` | `te_su`/`te_pr` aft corners disagree | station TE nodes |
| `mid_surface_fold` | Large skin–skin dihedral | `edges.csv` dihedral_deg |
| `edge_over_target` | Edge longer than recipe sizing | `edges.csv` vs `recipe.edge` |

## Lookalikes

| Do not merge | Why |
|---|---|
| `ShearWeb` / `WebWall` / `ShearWebComponent` | Station DTO; section recipe; SDF IR |
| `CellPolygon` / `CellWallLoop` | Loft cut interior; station FE loop DTO |
| `SectionWall` / `ShellMidlineStrip` | Solid midline on the cut; gated emit strip in \(\mathrm{B}\) |
| `cut.walls` solids / `shell_graph` partitions | Caps listed as solids; conformal skin-mid loop with cap colour |
| `thickness_at` / `query()["thickness"]` / `StackThickness` | Skin+ply field; stack \(T\) for `cavity_sdf` — do not merge |
| `build_shell_midline_strips` / `build_solid_cap_strips` / `shell_strips_from_graph` | Station skin+web; solid caps; conformal graph strips |
| `CapPlacementRule` / `AnchoredWallGeometry` | Plan-form cap policy; through-thickness face offsets |
| `CapPlacementKind` | Kind on `CapPlacementRule`, not a graph alias |
| `SectionCut` / `SectionGeometryResult` | Loft wrap at \(z\); gated station worker product |
| `SectionCut.layout` / `station_webs` / `station_family` | Authored concept token; realised W and X after web drop |

## Barrel functions

| Name | Role | Home |
|---|---|---|
| `loft_mapped_strips` | Map authored regions along the span | [blade](blade/README.md) |
| `rib_face_polygons` | Rib n-gons at a transverse station | [blade](blade/README.md) |
| `station_network` | Build `StationNetwork` | [blade](blade/README.md) |
| `wall_meta_from_cut` | Feature meta from a cut | [blade](blade/README.md) |
| `read_planform_csv` / `write_windio` | Planform I/O | [blade](blade/README.md) |
| `transverse_members_from_windio` | Ribs from a windIO structure | [blade](blade/README.md) |

Demoted from the barrel in 0.9.x (still on [blade](blade/README.md)):
`partition_edge_curve`, `shared_loft_names`, `signatures_compatible`,
`shell_strips_from_graph`, `station_network_from_meta`, `zone_region_graphs`,
`topology_diff`, `refine_section_samples`. Station helpers
`build_shell_midline_strips` / `build_solid_cap_strips` /
`emit_constructive_polygons` stay on their packages, not the barrel.
