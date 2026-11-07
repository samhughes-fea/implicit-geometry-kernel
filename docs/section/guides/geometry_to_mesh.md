# Geometry to mesh

Turn gated midlines into `section_mesh` inputs. Overview:
[`../README.md`](../docs_index.md). Runnable twin:
`../examples/05_handoff_to_mesh.py`.

↑ [`section_geometry`](../docs_index.md) · [`docs index`](../docs_index.md)

## Table of contents

1. [Pipeline](#pipeline)
2. [Required geometry gates](#required-geometry-gates)
3. [Adapter](#adapter)
4. [What `ShellMeshInputs` carries](#what-shellmeshinputs-carries)
5. [Optional coarse strip mesh](#optional-coarse-strip-mesh)
6. [Oracle policy](#oracle-policy-wave-0)
7. [What next](#what-next)

---

## Pipeline

```text
section (midlines) → ShellMidlineStrip list / shell_graph()
                   → downstream mesher
```

This package stays material-bind agnostic: handoff is wall centrelines and
thicknesses only.

## Spanwise geometry for a shell mesher (2.5-D)

Isolated `section_at(z)` midlines remain the layout oracle. Solids stay split
(`cut.walls` still lists caps). The **shell graph** is one wall: consume
`mid_surface()` / `shell_graph()` plus the live fields, not a second cap surface.
Sandwich cores and trailing-edge inserts stay constructive solids (not loop runs).

| Kernel | Role |
|---|---|
| `cut.shell_graph()` | Loftable skin loop (caps as partitions) + web branches; combined `T` |
| `graph.partition_edges()` / `partition_edge_curve` | Named feature lines `(z, name_a, name_b) → xy`; omit if a name is missing |
| `topology_diff` / `shared_loft_names` | Cap runout, dropped web, bay split; loft intersection only (no lerp) |
| `shell_strips_from_graph(cut)` | Shell `ShellMidlineStrip`s from the graph (no extra cap loop) |
| `ShellGraph.with_web_feet_on_loop()` | Web endpoints on the stacked loop (applied by `shell_graph`) |
| `BladeGeometry.section_at(z)` | Cut at `z`; LRU 256 keyed by `round(z, 9)` (no BVH) |
| `web_positions_from_windio_structure` | Chord fractions from windIO `structure.anchors` `nd_arc` |
| `geo.oml_sdf()` | Exact polygon SDF at each point's own `z` |
| `geo.cavity_sdf()` / `geo.wall_mid_sdf()` | cavity is `oml + T`; mid is the signed SDF of the `shell_graph` skin loop |
| `geo.wall_meta()` | Tagged bands of that graph; `evaluate((x,y,z)).meta` is the partition id |
| `GeometrySolver.project(..., plane=True)` | Snap nodes onto `phi = 0` in the station plane |
| `station_network(cut)` | Discrete graph from `shell_graph()` (loop runs + webs) |
| `station_network_from_meta(wall_mid, wall_meta, z, bounds, table)` | Samples the stacked mid loop, not a nested cap |
| `extract.feature_curve(a, b, zs, guess)` | Spanwise junction locator (`phi_a = phi_b = 0`) |

Do **not** copy Gmsh ThruSections, element `feature` arrays, or ply stacks into this
kernel. Colouring and conformal loft belong downstream, querying these fields.

---

## Oracle policy (Wave 0)

Two oracles, on purpose:

| Use | Source |
|---|---|
| Loft / conformal strips | Discrete `cut.shell_graph()` polylines |
| In-plane snap | `geo.wall_mid_sdf()` + `GeometrySolver.project(..., plane=True)` (`snap_wall_mid`) |

`T` is a **step**: skin gauge everywhere, `t_skin + t_cap` inside cap footprints
(`StackThickness`). No geometric taper (not ply stagger). On the 2D-F inventory
fixture that jump is **8 mm**; a shared cap/skin vertex may evaluate as either
3 mm or 11 mm.

`wall_mid_sdf` is the signed polygon SDF of that skin loop, so partition vertices
are on `φ = 0` as well as run interiors. `StackThickness` is still a **step**
(8 mm cap jump on the 2D-F fixture); that jump is a thickness oracle, not the
mid-surface field.

Numbers and remaining waves: [`MERGED_SHELL_GAPS.md`](../reports/MERGED_SHELL_GAPS.md).

`ShellMidlineStrip` is unchanged. Shell loft uses `shell_strips_from_graph` /
`shell_graph` (no extra cap surface). If a run name is missing at a station,
`shared_loft_names` omits that strip (`topology_diff` records the event).

---

## Required geometry gates

Include `"midlines"` (or `"cell_wall_loops"`, which also fills midlines):

```python
result = SectionGeometry(
    layout=resolve_system_type("2D-F"),
    airfoil=airfoil,
    skin_thickness=0.004,
    spar_cap_thickness=0.012,
    shear_web_thickness=0.005,
    outputs=("midlines",),
).execute().get_results()
```

Layouts `"0A"` / `"0B"` are valid shell-mesh stations (skin ± caps; no web
junctions). See
[`ADR_ZERO_WEB_SHELL_MESH.md`](../reference/ADR_ZERO_WEB_SHELL_MESH.md) and
[`feasibility_and_drops.md`](feasibility_and_drops.md#zero-web-keys).

---

## Adapter

Shell consumers in this kernel:

```python
from implicit_geometry_kernel import shell_strips_from_graph

strips = shell_strips_from_graph(cut)  # one strip per ShellRun; no extra cap surface
```

`cut.walls` still lists solid caps. `ShellMidlineStrip` is unchanged.

The out-of-repo solids path (`build_shell_midline_strips` →
`shell_mesh_inputs_from_geometry_result`) still emits separate cap strips until
that adapter is pointed here.

---

## What `ShellMeshInputs` carries

| Field | Role |
|-------|------|
| `midlines` | `ShellMidlineStrip` list. Shell path: `shell_strips_from_graph` (partitions). Solids path: still includes separate cap strips. |
| `chord_m`, `twist_rad` | Station scalars |
| `leading_edge_…`, `trailing_edge_…` | Twisted LE/TE markers in the blade frame |
| `layout_key` | Optional catalogue key for downstream logging |

Strip `kind` / `label` match the anatomy table in
[`anatomy_on_a_figure.md`](anatomy_on_a_figure.md)
(same folder as this guide).

---

## Optional coarse strip mesh

The kernel does not build a strip mesh. Consume `ShellMidlineStrip` (and, on a
loft, `shell_graph()`) in the downstream mesher.

Example flag: `python .../05_handoff_to_mesh.py --mesh`.

---

## What next

| Goal | Go to |
|------|--------|
| Midline emit details | [`../midlines/`](../midlines/README.md) |

Waves 0–7 are this page. G10/G11 stay parked in
[`MERGED_SHELL_GAPS.md#waves`](../reports/MERGED_SHELL_GAPS.md#waves).
Do not treat cores or inserts as missing from `shell_graph()`.
