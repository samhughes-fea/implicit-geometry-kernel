# Glossary

↑ [docs](README.md)

Look-up of public types. One home page per name. Barrel names are those in
`implicit_geometry_kernel.__all__`. Station and contract types are imported from
their packages.

Do not merge lookalike names.

## Contents

1. [Loft](#loft)
2. [Station](#station)
3. [Contract](#contract)
4. [SDF](#sdf)
5. [Lookalikes](#lookalikes)
6. [Barrel functions](#barrel-functions)

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
| `SectionWall` | One wall midline in the blade frame | `blade/cut.py` |
| `CellPolygon` | Closed cell interior on a loft cut | `blade/cut.py` |
| `ShellGraph` | Loftable skin loop, web branches, rib faces | `blade/shell_graph.py` |
| `ShellRun` | One loftable run on the graph | `blade/shell_graph.py` |
| `MidSurface` | Oriented mid-surface segments of a cut | `blade/shell_graph.py` |
| `PartitionEdge` | Meeting vertex between two loop partitions | `blade/shell_graph.py` |
| `shell_graph` | Builder: `SectionCut` → `ShellGraph` | `blade/shell_graph.py` |
| `BladeConcept` | Catalogue or parametric layout plus spanwise extras | `blade/concepts/concept.py` |
| `SpanwiseScheme` | Ribs, stringers, ply drops compiled to kwargs | `blade/concepts/concept.py` |
| `LockedLoftPlan` | Held point count per run across topology events | `blade/loft/locked_loft.py` |
| `resample_locked` | Resample graph runs to a locked plan | `blade/loft/locked_loft.py` |
| `shell_loft_zones` | Name-signature bands along the span | `blade/loft/locked_loft.py` |
| `WallRegionGraph` | Authored loop regions and webs | `blade/topology/graph.py` |
| `Region` | Named loop partition on the topology graph | `blade/topology/graph.py` |
| `Anchor` | Named vertex on a region graph | `blade/topology/graph.py` |
| `WebMember` | Web on the region graph | `blade/topology/graph.py` |
| `ShellTopologyTemplate` | Authored topology template | `blade/topology/template.py` |
| `StationNetwork` | Merged wall graph from feature meta | `blade/topology/network.py` |
| `TopologyEvent` | Onset or drop of a named run | `blade/shell_graph.py` |
| `TopologyDiff` | Difference of two graph signatures | `blade/shell_graph.py` |

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
| `LoftedPolygonSDF` | Polygon SDF at each point's `z`; `accelerate` BVH via `accel_at` / `bvh` |
| `GeometrySolver` | In-plane Newton (`sdf.solve`); not the extract package |

## Lookalikes

| Do not merge | Why |
|---|---|
| `ShearWeb` / `WebWall` / `ShearWebComponent` | Station DTO; section recipe; SDF IR |
| `CellPolygon` / `CellWallLoop` | Loft cut interior; station FE loop DTO |
| `SectionWall` / `ShellMidlineStrip` | Solid midline on the cut; gated emit strip in \(\mathrm{B}\) |
| `cut.walls` solids / `shell_graph` partitions | Caps listed as solids; conformal skin-mid loop with cap colour |
| `CapPlacementRule` / `AnchoredWallGeometry` | Plan-form cap policy; through-thickness face offsets |
| `CapPlacementKind` | Kind on `CapPlacementRule`, not a graph alias |
| `SectionCut` / `SectionGeometryResult` | Loft wrap at \(z\); gated station worker product |
| `SectionCut.layout` / `station_webs` / `station_family` | Authored concept token; realised W and X after web drop |

## Barrel functions

| Name | Role | Home |
|---|---|---|
| `loft_mapped_strips` | Map authored regions along the span | [blade](blade/README.md) |
| `partition_edge_curve` | Span curve of a partition edge | [blade](blade/README.md) |
| `rib_face_polygons` | Rib n-gons at a transverse station | [blade](blade/README.md) |
| `shared_loft_names` | Names common to two graphs | [blade](blade/README.md) |
| `signatures_compatible` | Whether two topology signatures loft | [blade](blade/README.md) |
| `shell_strips_from_graph` | Graph runs → contract strips | [geometry to mesh](section/guides/geometry_to_mesh.md) |
| `station_network` / `station_network_from_meta` | Build `StationNetwork` | [blade](blade/README.md) |
| `wall_meta_from_cut` | Feature meta from a cut | [blade](blade/README.md) |
| `zone_region_graphs` | Region graphs per loft zone | [blade](blade/README.md) |
| `topology_diff` | Diff two graph signatures | [blade](blade/README.md) |
| `emit_constructive_polygons` | Wall solids + sandwich cores as polygons (hybrid anatomy / example 09 `csg_mesh`) | `section/boundary/cores.py` |
| `read_planform_csv` / `write_windio` | Planform I/O | [blade](blade/README.md) |
| `transverse_members_from_windio` | Ribs from a windIO structure | [blade](blade/README.md) |
