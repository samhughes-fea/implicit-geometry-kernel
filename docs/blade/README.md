# blade

Span loft: airfoils, stations, `BladeGeometry`, cuts, and the merged shell graph.

**Import:** `from implicit_geometry_kernel import BladeGeometry, Station, Airfoil`
(also `implicit_geometry_kernel.blade`).

↑ [docs](../README.md) · [src](../../src/implicit_geometry_kernel/blade/README.md) ·
[examples](../../examples/blade/README.md) · [tests](../../tests/blade/README.md)

## What it can do

- Define stations (`Station`: span, chord, twist, pitch axis, offset, `Airfoil`).
- Loft a layout (`BladeGeometry`) and cut anywhere (`section_at` → `SectionCut`).
- Query OML and wall oracles: `oml_polyline`, `oml_sdf`, `wall_mid_sdf`, `wall_meta`.
- Merge wall midlines into `shell_graph` / `mid_surface()` for loftable strips.
- Read planform / windIO helpers; write OML STL and mid-surface patches.

`section_at` calls `SectionGeometry`. Loft oracles are extra fields, not a second
geometry engine.

## Mapped loft freeze

When `dropped_web_indices` is passed, loft webs equal intercept
`DroppedWeb.index` (no 0.08 cluster, no extra TE-plug omit).
`_regularise_anchors` keeps `web_*` over a nearby `cap_*`. TE wrap includes the
OML ridge; TE metrics use ridge / max chord-frac, not max blade-x. Placement
graph law equals the section law; skin loft at a web event is the t0 mould walk.
Do not add per-blade loft branches. Legacy `_cluster_web_members` runs only when
`dropped_web_indices is None`. Standing rules and closed loft bugs:
[bug_log.md](bug_log.md). Consumer table:
[geometry_to_mesh.md](../section/guides/geometry_to_mesh.md#what-next).

## Public names

`Airfoil`, `Station`, `SkinRegion`, `ShearWeb`, `Stringer`, `TransverseMember`,
`PlyDrop`, `BladeGeometry`, `BladeConcept`, `SectionCut`, `SectionWall`,
`CellPolygon`, `shell_graph`, `ShellGraph`, `shell_loft_zones`, `LockedLoftPlan`,
`resample_locked`, `write_oml` / `write_mid_surfaces` (on `BladeGeometry`).

## Examples

| Folder | Shows |
|---|---|
| [01_loft_and_cut](../../examples/blade/01_loft_and_cut/) | Tapered twisted loft, cuts |
| [02_web_placement](../../examples/blade/02_web_placement/) | Chord fractions vs offsets; dropped webs |
| [03_airfoils_and_root_transition](../../examples/blade/03_airfoils_and_root_transition/) | NACA, ellipse, flatback; circular root |
| [04_airfoil_families](../../examples/blade/04_airfoil_families/) | CST, t/c family, NACA 5/6 |
| [05_pchip_span](../../examples/blade/05_pchip_span/) | Monotone cubic planform |
| [06_blade_queries](../../examples/blade/06_blade_queries/) | Thickness field, batched query, validate |
| [07_mid_surfaces](../../examples/blade/07_mid_surfaces/) | Lofted mid-surface OBJ next to OML |
| [08_reference_blades](../../examples/blade/08_reference_blades/) | windIO 2.0 blades |
| [09_blade_concepts](../../examples/blade/09_blade_concepts/) | Catalogue concepts on a synthetic planform |

## See also

- One-station worker: [section](../section/README.md)
- Shell hand-off: [contract](../contract/README.md)
- Shell-graph gaps: [section/reports/MERGED_SHELL_GAPS.md](../section/reports/MERGED_SHELL_GAPS.md)
- Loft bug log: [bug_log.md](bug_log.md)
- Geometry-to-mesh handoff: [geometry_to_mesh.md](../section/guides/geometry_to_mesh.md)
