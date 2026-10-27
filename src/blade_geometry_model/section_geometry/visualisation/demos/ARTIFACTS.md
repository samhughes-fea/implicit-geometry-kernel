# Demo artefacts (regen index)

Figures under `../output/` (gitignored). CDFAM deck copies live in
`docs/conferences/CDFAM-2026/presentation/images/geometry/` via
`docs/conferences/CDFAM-2026/presentation/sync-assets.py`.

Run demos from the repository root. Prefix below is
`python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/`.

| Artefact | Demo command |
|----------|----------------|
| `bug_log_sg009/bug_log_sg009_schematic.png` | `.../bug_log_sg009_demo.py --presentation` |
| `feasibility_cascade/feasibility_cascade.png` | `.../feasibility_cascade_demo.py --presentation` |
| `csg_tree_diff/csg_tree_diff.png` | `.../csg_tree_diff_demo.py --presentation` |
| `structural_family_grids/structural_family_wall.png` | `.../structural_family_grids_demo.py --presentation --wall` |
| `anatomy_2d_f/2D-F_anatomy_annotated.png` | `.../anatomy_2d_f_demo.py` |
| `brep_vs_implicit_offset/brep_vs_implicit_offset.png` | `.../brep_vs_implicit_offset_demo.py --presentation` |
| `design_vector_morph/design_vector_morph.gif` | `.../design_vector_morph.py --presentation` |
| `design_vector_morph/design_vector_morph_stages.png` | `.../design_vector_morph.py --presentation` |
| `structural_family_grids/0A.png` | `.../structural_family_grids_demo.py` |
| `structural_family_grids/2D-F.png` | `.../structural_family_grids_demo.py` |
| `structural_family_grids/3B-F.png` | `.../structural_family_grids_demo.py` |
| `structural_family_grids/5A-F.png` | `.../structural_family_grids_demo.py` |
| `csg_expression_tree/csg_expression_tree.png` | `.../csg_expression_tree_demo.py` |
| `csg_expression_tree/csg_expression_tree_stages.gif` | `.../csg_expression_tree_demo.py --presentation` |
| `csg_expression_tree/csg_expression_tree_stages_trio.png` | `.../csg_expression_tree_demo.py --presentation` |
| `csg_expression_tree/csg_node_census.png` | `.../csg_expression_tree_demo.py` |
| `sdf_operator_algebra/sdf_operator_algebra.png` | `.../sdf_operator_algebra_demo.py` |
| `sdf_operator_algebra/sdf_kernel_only_operators.png` | `.../sdf_operator_algebra_demo.py` |
| `thickness_as_field_offset/thickness_sweep.png` | `.../thickness_as_field_offset_demo.py` |
| `thickness_as_field_offset/offset_metric_validity.png` | `.../thickness_as_field_offset_demo.py` |
| `sdf_metric_quality/metric_quality_fields.png` | `.../sdf_metric_quality_demo.py` |
| `sensitivity_dashboard/sensitivity_slide.png` | `.../sensitivity_dashboard_demo.py --presentation` |
| `le_te_exactness/le_nose_grid_vs_exact.png` | `.../le_te_exactness_demo.py` |
| `le_trailing_edge_insert_residual/inserts_le_te_zooms.png` | `.../le_te_insert_residual_demo.py --presentation` |
| `csg_provenance/csg_provenance_map.png` | `.../csg_provenance_demo.py` |
| `dual_rep_consistency/dual_rep_boundary_overlay.png` | `.../dual_rep_consistency_demo.py` |
| `dual_rep_consistency/dual_rep_residual_track.png` | `.../dual_rep_consistency_demo.py` |
| `thin_feature_robustness/thin_web_zoom_twist.png` | `.../thin_feature_robustness_demo.py` |
| `robustness_monte_carlo/robustness_monte_carlo.png` | `.../robustness_monte_carlo_demo.py` |
| `quadtree_union/quadtree_union_5C-F.png` | `.../quadtree_union_mesh.py` |
| `evaluation_cost_pruning/support_box_culling.png` | `.../evaluation_cost_pruning_demo.py` |
| `evaluation_cost_pruning/evaluation_scaling.png` | `.../evaluation_cost_pruning_demo.py` |
| `props_cost/props_cost_pareto.png` | `.../props_cost_demo.py` |
| `manystation_bench/manystation_throughput.png` | `.../manystation_bench.py --stations 12 --figure --presentation` |
| `gifs/gif_02_system_taxonomy_40keys.gif` | `.../feature_flythroughs.py` |
| `spanwise_yz_plan_flyover/spanwise_yz_plan_flyover.gif` | `.../spanwise_yz_plan_flyover.py` |
| `topology_twist_stability/gap_vs_twist.png` | `.../topology_twist_stability_demo.py` |
| `shell_junction_topology/topology_junction_zooms.png` | `.../shell_junction_topology_demo.py` |
| `geometry_to_mesh/geometry_to_mesh_flythrough.gif` | `.../geometry_to_mesh_flythrough.py` |
| `geometry_to_mesh/geometry_to_mesh_stages.png` | `.../geometry_to_mesh_flythrough.py` |
| `cell_wall_loops/cell_wall_loops_on_section_n_webs_2.png` | `.../cell_wall_loops_demo.py --presentation` |
| `api_quickstart/api_quickstart_section.png` | `.../api_quickstart_demo.py --presentation` |

**Built by `sync-assets.py` (not demos):** after copying the rows above into
`docs/conferences/CDFAM-2026/presentation/images/geometry/`, the script also
writes:

| Deck image | How |
|------------|-----|
| `taxonomy_still.png` | 2×2 tile of `0A.png`, `2D-F.png`, `3B-F.png`, `5A-F.png` |
| `flyover_still.png` | frame ~80% through `spanwise_yz_plan_flyover.gif` |

Chrome / profile / WP images are copied from the ICOE deck (see `_from_icoe` in
`sync-assets.py`), not from `section_geometry` demos.
