# visualisation demos — run list

Figure gallery for station geometry. This is not the first on-ramp: start with
[`../../section/teaching/`](../../section/teaching/) then the numbered examples
in [`../../README.md`](../../README.md). Outputs:
[`../figures/`](../figures/), committed and shown under [Figures](#figures) below. Artefact / regen index:
[`ARTIFACTS.md`](ARTIFACTS.md). Run from the repository root.

See also [`../README.md`](../README.md).

```bash
python examples/demos/section/<script>.py
```

Scripts listed in [`scripts/run_examples.py`](../../../scripts/run_examples.py)
`SKIP` stay in-tree but are not run in CI:
[`le_te_insert_residual_demo.py`](le_te_insert_residual_demo.py) (hybrid residual
quadtree no longer exposed) and [`quadtree_union_mesh.py`](quadtree_union_mesh.py)
(full-resolution, several minutes).

## Table of contents

1. [Gallery tour](#gallery-tour)
2. [Product](#product)
3. [Exactness](#exactness)
4. [Kernel](#kernel)
5. [Perf / benches](#perf--benches)
6. [Stdout-only / helpers](#stdout-only--helpers)
7. [Not ported yet](#not-ported-yet)
8. [Figures](#figures)

---

## Gallery tour

After teaching `01`–`06`, these gallery scripts are a useful figure tour
(from the repository root):

1. [`naca2412_section.py`](naca2412_section.py) — one twin-web station
2. [`structural_family_grids_demo.py`](structural_family_grids_demo.py) — all 40 layout keys
3. [`anatomy_2d_f_demo.py`](anatomy_2d_f_demo.py) — labelled `2D-F` anatomy
4. [`shell_junction_topology_demo.py`](shell_junction_topology_demo.py) — junction snap
5. [`cell_wall_loops_demo.py`](cell_wall_loops_demo.py) — cell-wall loops

Checked-in gallery:
[`../../../docs/section/reference/structural_systems_gallery.md`](../../../docs/section/reference/structural_systems_gallery.md).
Regen commands: [`ARTIFACTS.md`](ARTIFACTS.md).

Figures use [`figure_io.py`](../../../src/implicit_geometry_kernel/plotting/figure_io.py)
(`save_figure`, default DPI 150). Conference-facing demos accept `--presentation`
(larger fonts, DPI 220).

---

## Product

Section / gallery / junctions / cell walls / anatomy.

| Script | What it draws |
|--------|----------------|
| [`naca2412_section.py`](naca2412_section.py) | Twin-web section and twisted box examples — `../figures/naca2412_section/` |
| [`geometry_from_design_vector.py`](geometry_from_design_vector.py) | Build/plot from thicknesses + layout (families A-D) |
| [`structural_family_grids_demo.py`](structural_family_grids_demo.py) | One PNG per layout id; `--wall` contact sheet — `../figures/structural_family_grids/` |
| [`anatomy_2d_f_demo.py`](anatomy_2d_f_demo.py) | Labelled `2D-F` with legend — `../figures/anatomy_2d_f/` |
| [`section_centerline_geometry_taxonomy.py`](section_centerline_geometry_taxonomy.py) | Per-label SDF panels for layout keys |
| [`system_types_twist_grid.py`](system_types_twist_grid.py) | Layout grids at fixed twist |
| [`skin_ss_ps_taxonomy_demo.py`](skin_ss_ps_taxonomy_demo.py) | Upper/lower skin solids and midline arcs |
| [`shell_junction_topology_demo.py`](shell_junction_topology_demo.py) | Junction snap before/after (raw vs prepared midlines) |
| [`cell_wall_loops_demo.py`](cell_wall_loops_demo.py) | Closed cell loops for 0-5 webs |
| [`api_quickstart_demo.py`](api_quickstart_demo.py) | The README quickstart, run, plus what the result holds |

---

## Exactness

| Script | What it draws |
|--------|----------------|
| [`sg009_axis_cap_junction_spike.py`](sg009_axis_cap_junction_spike.py) | SG-009 offline: axis ? cap oracle vs production shared nodes |
| [`csg_tree_diff_demo.py`](csg_tree_diff_demo.py) | Diff two CSG trees for one region label |
| [`feasibility_cascade_demo.py`](feasibility_cascade_demo.py) | Chord-shrink resolve ladder: drop / step-down / raise |
| [`bug_log_sg009_demo.py`](bug_log_sg009_demo.py) | SG-009 tangency schematic (crossing then shared node) |
| [`dual_rep_consistency_demo.py`](dual_rep_consistency_demo.py) | Polygon boundary vs SDF contour |
| [`le_te_exactness_demo.py`](le_te_exactness_demo.py) | Skin midline vs aerofoil tip |
| [`le_te_insert_residual_demo.py`](le_te_insert_residual_demo.py) | LE/TE inserts on the SDF path (in-tree; **CI SKIP**) |
| [`sdf_metric_quality_demo.py`](sdf_metric_quality_demo.py) | SDF gradient quality fields |
| [`thin_feature_robustness_demo.py`](thin_feature_robustness_demo.py) | Thin webs under twist |
| [`topology_twist_stability_demo.py`](topology_twist_stability_demo.py) | Junction gaps vs twist |
| [`robustness_monte_carlo_demo.py`](robustness_monte_carlo_demo.py) | Randomized thickness+twist robustness sweep |
| [`sensitivity_dashboard_demo.py`](sensitivity_dashboard_demo.py) | Finite-difference sensitivities and local sweeps |
| [`csg_provenance_demo.py`](csg_provenance_demo.py) | Winner-label / nearest-interface / ambiguity maps |
| [`brep_vs_implicit_offset_demo.py`](brep_vs_implicit_offset_demo.py) | Curve offset folds; the level set cannot |

---

## Kernel

Implicit-geometry machinery (presentation order: expression tree ? operator
algebra ? thickness-as-offset ? evaluation cost ? geometry-to-mesh fly-through).

Stages 2 and 3 are both "midlines" and are commonly conflated: **raw** emit has
free endpoints; **prepared** topology splits so junctions share nodes. Meshing
consumes the prepared form only.

| Script | What it draws |
|--------|----------------|
| [`csg_expression_tree_demo.py`](csg_expression_tree_demo.py) | The `Expr` tree behind one region |
| [`sdf_operator_algebra_demo.py`](sdf_operator_algebra_demo.py) | Every IR operator as one line of numpy |
| [`thickness_as_field_offset_demo.py`](thickness_as_field_offset_demo.py) | Gauge sweep / Eikonal residual |
| [`evaluation_cost_pruning_demo.py`](evaluation_cost_pruning_demo.py) | Shared-subtree memoisation and bbox cull |

---

## Perf / benches

| Script | What it draws |
|--------|----------------|
| [`props_cost_demo.py`](props_cost_demo.py) | Hybrid vs quadtree vs grid cost |
| [`quadtree_union_mesh.py`](quadtree_union_mesh.py) | Shared union quadtree meshes (in-tree; **CI SKIP**, slow) |
| [`quadtree_convergence_study.py`](quadtree_convergence_study.py) | Quadtree depth vs dense-grid reference |
| [`manystation_bench.py`](manystation_bench.py) | Throughput across many stations |

---

## Stdout-only / helpers

| Script | Role |
|--------|------|
| [`constructive_vs_quadtree_spike.py`](constructive_vs_quadtree_spike.py) | Console comparison (stdout) |
| [`_ensure_repo_path.py`](_ensure_repo_path.py) | Put repo root on `sys.path` |
| [`_taxonomy_cases.py`](_taxonomy_cases.py) | Shared twist/aerofoil helpers |
| [`cell_wall_loops_case_payload.py`](cell_wall_loops_case_payload.py) | JSON payload helper |
| [`_plot_debug.py`](_plot_debug.py) | Optional Matplotlib call tracing |

Mesh strip demos belong in a downstream mesher, not this kernel.

Dual-rep acceptance narrative:
[`../../../docs/section/reports/DUAL_REP_ACCEPTANCE.md`](../../../docs/section/reports/DUAL_REP_ACCEPTANCE.md).

---

## Not ported yet

Fork-era drivers cited in older docs. Names kept here so the gallery index is
honest; do not look for these `.py` files in this folder.

| Name | Role (when it existed) |
|------|------------------------|
| `spanwise_web_anchors.py` | Root-anchored vs per-station web positions (see teaching `04`) |
| `section_properties_dashboard.py` | Per-label area / inertia dashboard |
| `feature_flythroughs.py` | Reusable GIF fly-through pipeline |
| `design_vector_morph.py` | Design-vector morph GIF (static plot: [`geometry_from_design_vector.py`](geometry_from_design_vector.py)) |
| `geometry_to_mesh_flythrough.py` | Solids → midlines → Q4 mesh (kernel stops at contract types) |
| `spanwise_yz_plan_flyover.py` | Spanwise YZ plan fly-over GIF |
| `quadtree_vs_shoelace_convergence_study.py` | Polygon vs quadtree accuracy–cost (see [`../../integration/01_quadtree_vs_shoelace/`](../../integration/01_quadtree_vs_shoelace/)) |

### Plot debug tracing

Set `SECTION_GEOMETRY_PLOT_DEBUG=1` to capture structured summaries of data
passed to Matplotlib calls for any demo script.

- Default log dir: `../figures/_plot_debug/`
- Optional overrides: `SECTION_GEOMETRY_PLOT_DEBUG_DIR`, `SECTION_GEOMETRY_PLOT_DEBUG_TAG`

---

## Figures

Every figure the demos write, committed under [`../figures/`](../figures/) so it renders here. `python scripts/run_examples.py` regenerates them.

### `naca2412_section` · [`naca2412_section.py`](naca2412_section.py)

![component_sdfs_torsion](../figures/naca2412_section/component_sdfs_torsion.png)

![torsion_box_section](../figures/naca2412_section/torsion_box_section.png)

![torsion_box_section_chord_normal](../figures/naca2412_section/torsion_box_section_chord_normal.png)

![torsion_box_section_flapwise](../figures/naca2412_section/torsion_box_section_flapwise.png)

![twin_web_section](../figures/naca2412_section/twin_web_section.png)

### `structural_family_grids` · [`structural_family_grids_demo.py`](structural_family_grids_demo.py)

<table>
<tr><td align="center"><img src="../figures/structural_family_grids/0A.png" width="200"><br><sub>0A</sub></td><td align="center"><img src="../figures/structural_family_grids/0B.png" width="200"><br><sub>0B</sub></td><td align="center"><img src="../figures/structural_family_grids/1A-CN.png" width="200"><br><sub>1A-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/1A-F.png" width="200"><br><sub>1A-F</sub></td></tr>
<tr><td align="center"><img src="../figures/structural_family_grids/1B-CN.png" width="200"><br><sub>1B-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/1B-F.png" width="200"><br><sub>1B-F</sub></td><td align="center"><img src="../figures/structural_family_grids/1C-CN.png" width="200"><br><sub>1C-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/1C-F.png" width="200"><br><sub>1C-F</sub></td></tr>
<tr><td align="center"><img src="../figures/structural_family_grids/2A-CN.png" width="200"><br><sub>2A-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/2A-F.png" width="200"><br><sub>2A-F</sub></td><td align="center"><img src="../figures/structural_family_grids/2B-CN.png" width="200"><br><sub>2B-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/2B-F.png" width="200"><br><sub>2B-F</sub></td></tr>
<tr><td align="center"><img src="../figures/structural_family_grids/2C-CN.png" width="200"><br><sub>2C-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/2C-F.png" width="200"><br><sub>2C-F</sub></td><td align="center"><img src="../figures/structural_family_grids/2D-CN.png" width="200"><br><sub>2D-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/2D-F.png" width="200"><br><sub>2D-F</sub></td></tr>
<tr><td align="center"><img src="../figures/structural_family_grids/3A-CN.png" width="200"><br><sub>3A-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/3A-F.png" width="200"><br><sub>3A-F</sub></td><td align="center"><img src="../figures/structural_family_grids/3B-CN.png" width="200"><br><sub>3B-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/3B-F.png" width="200"><br><sub>3B-F</sub></td></tr>
<tr><td align="center"><img src="../figures/structural_family_grids/3C-CN.png" width="200"><br><sub>3C-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/3C-F.png" width="200"><br><sub>3C-F</sub></td><td align="center"><img src="../figures/structural_family_grids/3D-CN.png" width="200"><br><sub>3D-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/3D-F.png" width="200"><br><sub>3D-F</sub></td></tr>
<tr><td align="center"><img src="../figures/structural_family_grids/4A-CN.png" width="200"><br><sub>4A-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/4A-F.png" width="200"><br><sub>4A-F</sub></td><td align="center"><img src="../figures/structural_family_grids/4B-CN.png" width="200"><br><sub>4B-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/4B-F.png" width="200"><br><sub>4B-F</sub></td></tr>
<tr><td align="center"><img src="../figures/structural_family_grids/4C-CN.png" width="200"><br><sub>4C-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/4C-F.png" width="200"><br><sub>4C-F</sub></td><td align="center"><img src="../figures/structural_family_grids/4D-CN.png" width="200"><br><sub>4D-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/4D-F.png" width="200"><br><sub>4D-F</sub></td></tr>
<tr><td align="center"><img src="../figures/structural_family_grids/5A-CN.png" width="200"><br><sub>5A-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/5A-F.png" width="200"><br><sub>5A-F</sub></td><td align="center"><img src="../figures/structural_family_grids/5B-CN.png" width="200"><br><sub>5B-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/5B-F.png" width="200"><br><sub>5B-F</sub></td></tr>
<tr><td align="center"><img src="../figures/structural_family_grids/5C-CN.png" width="200"><br><sub>5C-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/5C-F.png" width="200"><br><sub>5C-F</sub></td><td align="center"><img src="../figures/structural_family_grids/5D-CN.png" width="200"><br><sub>5D-CN</sub></td><td align="center"><img src="../figures/structural_family_grids/5D-F.png" width="200"><br><sub>5D-F</sub></td></tr>
</table>

### `anatomy_2d_f` · [`anatomy_2d_f_demo.py`](anatomy_2d_f_demo.py)

![2D-F_anatomy_annotated](../figures/anatomy_2d_f/2D-F_anatomy_annotated.png)

![2D-F_labelled](../figures/anatomy_2d_f/2D-F_labelled.png)

### `shell_junction_topology` · [`shell_junction_topology_demo.py`](shell_junction_topology_demo.py)

![topology_before_after](../figures/shell_junction_topology/topology_before_after.png)

![topology_junction_zooms](../figures/shell_junction_topology/topology_junction_zooms.png)

### `cell_wall_loops` · [`cell_wall_loops_demo.py`](cell_wall_loops_demo.py)

<table>
<tr><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_n_webs_0.png" width="200"><br><sub>cell_wall_loops_n_webs_0</sub></td><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_n_webs_1.png" width="200"><br><sub>cell_wall_loops_n_webs_1</sub></td><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_n_webs_2.png" width="200"><br><sub>cell_wall_loops_n_webs_2</sub></td><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_n_webs_3.png" width="200"><br><sub>cell_wall_loops_n_webs_3</sub></td></tr>
<tr><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_n_webs_4.png" width="200"><br><sub>cell_wall_loops_n_webs_4</sub></td><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_n_webs_5.png" width="200"><br><sub>cell_wall_loops_n_webs_5</sub></td><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_on_section_n_webs_0.png" width="200"><br><sub>cell_wall_loops_on_section_n_webs_0</sub></td><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_on_section_n_webs_1.png" width="200"><br><sub>cell_wall_loops_on_section_n_webs_1</sub></td></tr>
<tr><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_on_section_n_webs_2.png" width="200"><br><sub>cell_wall_loops_on_section_n_webs_2</sub></td><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_on_section_n_webs_3.png" width="200"><br><sub>cell_wall_loops_on_section_n_webs_3</sub></td><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_on_section_n_webs_4.png" width="200"><br><sub>cell_wall_loops_on_section_n_webs_4</sub></td><td align="center"><img src="../figures/cell_wall_loops/cell_wall_loops_on_section_n_webs_5.png" width="200"><br><sub>cell_wall_loops_on_section_n_webs_5</sub></td></tr>
</table>

### `design_vector` · [`geometry_from_design_vector.py`](geometry_from_design_vector.py)

![section_design_vector](../figures/design_vector/section_design_vector.png)

### `section_centerline_geometry_taxonomy` · [`section_centerline_geometry_taxonomy.py`](section_centerline_geometry_taxonomy.py)

<table>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_0A_i000.png" width="200"><br><sub>section_geom_0A_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_0B_i000.png" width="200"><br><sub>section_geom_0B_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_1A-CN_i000.png" width="200"><br><sub>section_geom_1A-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_1A-F_i000.png" width="200"><br><sub>section_geom_1A-F_i000</sub></td></tr>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_1B-CN_i000.png" width="200"><br><sub>section_geom_1B-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_1B-F_i000.png" width="200"><br><sub>section_geom_1B-F_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_1C-CN_i000.png" width="200"><br><sub>section_geom_1C-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_1C-F_i000.png" width="200"><br><sub>section_geom_1C-F_i000</sub></td></tr>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_2A-CN_i000.png" width="200"><br><sub>section_geom_2A-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_2A-F_i000.png" width="200"><br><sub>section_geom_2A-F_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_2B-CN_i000.png" width="200"><br><sub>section_geom_2B-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_2B-F_i000.png" width="200"><br><sub>section_geom_2B-F_i000</sub></td></tr>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_2C-CN_i000.png" width="200"><br><sub>section_geom_2C-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_2C-F_i000.png" width="200"><br><sub>section_geom_2C-F_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_2D-CN_i000.png" width="200"><br><sub>section_geom_2D-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_2D-F_i000.png" width="200"><br><sub>section_geom_2D-F_i000</sub></td></tr>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_3A-CN_i000.png" width="200"><br><sub>section_geom_3A-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_3A-F_i000.png" width="200"><br><sub>section_geom_3A-F_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_3B-CN_i000.png" width="200"><br><sub>section_geom_3B-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_3B-F_i000.png" width="200"><br><sub>section_geom_3B-F_i000</sub></td></tr>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_3C-CN_i000.png" width="200"><br><sub>section_geom_3C-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_3C-F_i000.png" width="200"><br><sub>section_geom_3C-F_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_3D-CN_i000.png" width="200"><br><sub>section_geom_3D-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_3D-F_i000.png" width="200"><br><sub>section_geom_3D-F_i000</sub></td></tr>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_4A-CN_i000.png" width="200"><br><sub>section_geom_4A-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_4A-F_i000.png" width="200"><br><sub>section_geom_4A-F_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_4B-CN_i000.png" width="200"><br><sub>section_geom_4B-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_4B-F_i000.png" width="200"><br><sub>section_geom_4B-F_i000</sub></td></tr>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_4C-CN_i000.png" width="200"><br><sub>section_geom_4C-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_4C-F_i000.png" width="200"><br><sub>section_geom_4C-F_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_4D-CN_i000.png" width="200"><br><sub>section_geom_4D-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_4D-F_i000.png" width="200"><br><sub>section_geom_4D-F_i000</sub></td></tr>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_5A-CN_i000.png" width="200"><br><sub>section_geom_5A-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_5A-F_i000.png" width="200"><br><sub>section_geom_5A-F_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_5B-CN_i000.png" width="200"><br><sub>section_geom_5B-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_5B-F_i000.png" width="200"><br><sub>section_geom_5B-F_i000</sub></td></tr>
<tr><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_5C-CN_i000.png" width="200"><br><sub>section_geom_5C-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_5C-F_i000.png" width="200"><br><sub>section_geom_5C-F_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_5D-CN_i000.png" width="200"><br><sub>section_geom_5D-CN_i000</sub></td><td align="center"><img src="../figures/section_centerline_geometry_taxonomy/section_geom_5D-F_i000.png" width="200"><br><sub>section_geom_5D-F_i000</sub></td></tr>
</table>

### `system_types_twist_grid` · [`system_types_twist_grid.py`](system_types_twist_grid.py)

![system_types_ABC_twist20deg](../figures/system_types_twist_grid/system_types_ABC_twist20deg.png)

### `skin_ss_ps_taxonomy` · [`skin_ss_ps_taxonomy_demo.py`](skin_ss_ps_taxonomy_demo.py)

![skin_ss_ps_midline_arcs](../figures/skin_ss_ps_taxonomy/skin_ss_ps_midline_arcs.png)

![skin_ss_ps_section](../figures/skin_ss_ps_taxonomy/skin_ss_ps_section.png)

### `api_quickstart` · [`api_quickstart_demo.py`](api_quickstart_demo.py)

![api_quickstart_section](../figures/api_quickstart/api_quickstart_section.png)

### `csg_tree_diff` · [`csg_tree_diff_demo.py`](csg_tree_diff_demo.py)

![csg_tree_diff](../figures/csg_tree_diff/csg_tree_diff.png)

### `feasibility_cascade` · [`feasibility_cascade_demo.py`](feasibility_cascade_demo.py)

![feasibility_cascade](../figures/feasibility_cascade/feasibility_cascade.png)

### `bug_log_sg009` · [`bug_log_sg009_demo.py`](bug_log_sg009_demo.py)

![bug_log_sg009_schematic](../figures/bug_log_sg009/bug_log_sg009_schematic.png)

### `dual_rep_consistency` · [`dual_rep_consistency_demo.py`](dual_rep_consistency_demo.py)

![dual_rep_boundary_overlay](../figures/dual_rep_consistency/dual_rep_boundary_overlay.png)

![dual_rep_gap_map](../figures/dual_rep_consistency/dual_rep_gap_map.png)

![dual_rep_residual_track](../figures/dual_rep_consistency/dual_rep_residual_track.png)

### `le_te_exactness` · [`le_te_exactness_demo.py`](le_te_exactness_demo.py)

![le_nose_grid_vs_exact](../figures/le_te_exactness/le_nose_grid_vs_exact.png)

![le_nose_sdf_fidelity](../figures/le_te_exactness/le_nose_sdf_fidelity.png)

![skin_midline_vs_aero_te](../figures/le_te_exactness/skin_midline_vs_aero_te.png)

![te_tip_grid_vs_hybrid](../figures/le_te_exactness/te_tip_grid_vs_hybrid.png)

![twist15_le_te_markers](../figures/le_te_exactness/twist15_le_te_markers.png)

### `sdf_metric_quality` · [`sdf_metric_quality_demo.py`](sdf_metric_quality_demo.py)

![eikonal_residual_hist](../figures/sdf_metric_quality/eikonal_residual_hist.png)

![metric_quality_fields](../figures/sdf_metric_quality/metric_quality_fields.png)

### `thin_feature_robustness` · [`thin_feature_robustness_demo.py`](thin_feature_robustness_demo.py)

![thin_web_zoom_twist](../figures/thin_feature_robustness/thin_web_zoom_twist.png)

![web_area_vs_depth](../figures/thin_feature_robustness/web_area_vs_depth.png)

### `topology_twist_stability` · [`topology_twist_stability_demo.py`](topology_twist_stability_demo.py)

![gap_vs_twist](../figures/topology_twist_stability/gap_vs_twist.png)

![junction_before_after](../figures/topology_twist_stability/junction_before_after.png)

![junction_zoom_worst_twist](../figures/topology_twist_stability/junction_zoom_worst_twist.png)

### `robustness_monte_carlo` · [`robustness_monte_carlo_demo.py`](robustness_monte_carlo_demo.py)

![robustness_monte_carlo](../figures/robustness_monte_carlo/robustness_monte_carlo.png)

### `sensitivity_dashboard` · [`sensitivity_dashboard_demo.py`](sensitivity_dashboard_demo.py)

![sensitivity_dashboard](../figures/sensitivity_dashboard/sensitivity_dashboard.png)

![sensitivity_slide](../figures/sensitivity_dashboard/sensitivity_slide.png)

### `csg_provenance` · [`csg_provenance_demo.py`](csg_provenance_demo.py)

![csg_provenance_map](../figures/csg_provenance/csg_provenance_map.png)

### `brep_vs_implicit_offset` · [`brep_vs_implicit_offset_demo.py`](brep_vs_implicit_offset_demo.py)

![brep_vs_implicit_offset](../figures/brep_vs_implicit_offset/brep_vs_implicit_offset.png)

### `csg_expression_tree` · [`csg_expression_tree_demo.py`](csg_expression_tree_demo.py)

![csg_expression_tree](../figures/csg_expression_tree/csg_expression_tree.png)

![csg_expression_tree_stages](../figures/csg_expression_tree/csg_expression_tree_stages.gif)

![csg_expression_tree_stages_trio](../figures/csg_expression_tree/csg_expression_tree_stages_trio.png)

![csg_node_census](../figures/csg_expression_tree/csg_node_census.png)

### `sdf_operator_algebra` · [`sdf_operator_algebra_demo.py`](sdf_operator_algebra_demo.py)

![sdf_kernel_only_operators](../figures/sdf_operator_algebra/sdf_kernel_only_operators.png)

![sdf_operator_algebra](../figures/sdf_operator_algebra/sdf_operator_algebra.png)

### `thickness_as_field_offset` · [`thickness_as_field_offset_demo.py`](thickness_as_field_offset_demo.py)

![offset_metric_validity](../figures/thickness_as_field_offset/offset_metric_validity.png)

![thickness_sweep](../figures/thickness_as_field_offset/thickness_sweep.png)

### `evaluation_cost_pruning` · [`evaluation_cost_pruning_demo.py`](evaluation_cost_pruning_demo.py)

![evaluation_scaling](../figures/evaluation_cost_pruning/evaluation_scaling.png)

![support_box_culling](../figures/evaluation_cost_pruning/support_box_culling.png)

### `props_cost` · [`props_cost_demo.py`](props_cost_demo.py)

![props_cost_pareto](../figures/props_cost/props_cost_pareto.png)

![props_cost_table](../figures/props_cost/props_cost_table.png)

### `quadtree_convergence` · [`quadtree_convergence_study.py`](quadtree_convergence_study.py)

![convergence_2web_t0](../figures/quadtree_convergence/convergence_2web_t0.png)

![convergence_2web_t15](../figures/quadtree_convergence/convergence_2web_t15.png)

### `brep_extend_trim` · [`brep_extend_trim_demo.py`](brep_extend_trim_demo.py)

![brep_extend_trim](../figures/brep_extend_trim/brep_extend_trim.png)

### `offset_shape_derivative` · [`offset_shape_derivative_demo.py`](offset_shape_derivative_demo.py)

![offset_shape_derivative](../figures/offset_shape_derivative/offset_shape_derivative.png)
