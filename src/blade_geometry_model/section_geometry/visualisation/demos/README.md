# visualisation demos ? run list

Demo scripts for `section_geometry` products. Outputs:
[`../output/`](../output/) (gitignored). Artefact / regen index:
[`ARTIFACTS.md`](ARTIFACTS.md). Run from the repository root.

Teaching scripts (print-first, few flags):
[`../../examples/`](../../examples/README.md).

? [`section_geometry`](../../README.md) ?
[`../README.md`](../README.md)

```bash
python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/<script>.py
```

## Table of contents

1. [Newcomer path](#newcomer-path)
2. [Product](#product)
3. [Exactness](#exactness)
4. [Kernel](#kernel)
5. [Perf / benches](#perf--benches)
6. [Stdout-only / helpers](#stdout-only--helpers)

---

## Newcomer path

Run these first (from the repository root):

1. [`naca2412_section.py`](naca2412_section.py) ? one twin-web station
2. [`structural_family_grids_demo.py`](structural_family_grids_demo.py) ? all 40 layout keys
3. [`anatomy_2d_f_demo.py`](anatomy_2d_f_demo.py) ? labelled `2D-F` anatomy
4. [`shell_junction_topology_demo.py`](shell_junction_topology_demo.py) ? junction snap
5. [`cell_wall_loops_demo.py`](cell_wall_loops_demo.py) ? cell-wall loops

Checked-in gallery:
[`../../docs/reference/structural_systems_gallery.md`](../../docs/reference/structural_systems_gallery.md).
Regen commands: [`ARTIFACTS.md`](ARTIFACTS.md).

Figures use [`../figure_io.py`](../figure_io.py) (`save_figure`, default DPI 200).
Conference-facing demos accept `--presentation` (larger fonts, DPI 220).

---

## Product

Section / gallery / junctions / cell walls / anatomy.

| Script | What it draws |
|--------|----------------|
| [`naca2412_section.py`](naca2412_section.py) | Twin-web section and twisted box examples ? `output/naca2412_section/` |
| [`geometry_from_design_vector.py`](geometry_from_design_vector.py) | Build/plot from thicknesses + layout (families A?D) |
| [`structural_family_grids_demo.py`](structural_family_grids_demo.py) | One PNG per layout id; `--wall` contact sheet ? `../output/structural_family_grids/` |
| [`anatomy_2d_f_demo.py`](anatomy_2d_f_demo.py) | Labelled `2D-F` with legend ? `../output/anatomy_2d_f/` |
| [`section_centerline_geometry_taxonomy.py`](section_centerline_geometry_taxonomy.py) | Per-label SDF panels for layout keys |
| [`system_types_twist_grid.py`](system_types_twist_grid.py) | Layout grids at fixed twist |
| [`spanwise_web_anchors.py`](spanwise_web_anchors.py) | Root-anchored vs per-station web positions |
| [`skin_ss_ps_taxonomy_demo.py`](skin_ss_ps_taxonomy_demo.py) | Upper/lower skin solids and midline arcs |
| [`shell_junction_topology_demo.py`](shell_junction_topology_demo.py) | Junction snap before/after (raw vs prepared midlines) |
| [`cell_wall_loops_demo.py`](cell_wall_loops_demo.py) | Closed cell loops for 0?5 webs |
| [`section_properties_dashboard.py`](section_properties_dashboard.py) | Per-label area / inertia dashboard |
| [`feature_flythroughs.py`](feature_flythroughs.py) | Reusable GIF fly-through pipeline |
| [`design_vector_morph.py`](design_vector_morph.py) | One design vector -> solids, midlines and Q4 mesh, animated |
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
| [`le_te_insert_residual_demo.py`](le_te_insert_residual_demo.py) | LE/TE inserts on the SDF path |
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
| [`geometry_to_mesh_flythrough.py`](geometry_to_mesh_flythrough.py) | Solids ? raw midlines ? prepared strips ? Q4 mesh |

The fly-through needs ImageMagick (`magick`) on PATH; pass `--no-gif` for
frames/static strip only. `--max-frames N` for a fast preview. `--presentation`
raises fonts/DPI.

---

## Perf / benches

| Script | What it draws |
|--------|----------------|
| [`props_cost_demo.py`](props_cost_demo.py) | Hybrid vs quadtree vs grid cost |
| [`quadtree_union_mesh.py`](quadtree_union_mesh.py) | Shared union quadtree meshes |
| [`quadtree_convergence_study.py`](quadtree_convergence_study.py) | Quadtree depth vs dense-grid reference |
| [`quadtree_vs_shoelace_convergence_study.py`](quadtree_vs_shoelace_convergence_study.py) | Polygon vs quadtree accuracy?cost |
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

Mesh strip demos live under
[`section_mesh/visualisation/demos/`](../../../../mesh_utilities/section_mesh/visualisation/demos/).

Dual-rep acceptance narrative:
[`../../docs/reports/DUAL_REP_ACCEPTANCE.md`](../../docs/reports/DUAL_REP_ACCEPTANCE.md).

### Plot debug tracing

Set `SECTION_GEOMETRY_PLOT_DEBUG=1` to capture structured summaries of data
passed to Matplotlib calls for any demo script.

- Default log dir: `../output/_plot_debug/`
- Optional overrides: `SECTION_GEOMETRY_PLOT_DEBUG_DIR`, `SECTION_GEOMETRY_PLOT_DEBUG_TAG`
