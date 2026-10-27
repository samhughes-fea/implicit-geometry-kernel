# Where this geometry came from

`section_geometry`, `section_contract`, `coordinate_frames`, `area_2d` and
`shoelace_integration` are a one-time fork, with history, of the same packages in
[samhughes-fea/blade-structure](https://github.com/samhughes-fea/blade-structure)
at commit `51ba12f0` (27 October 2026), where `section_geometry` is complete and
frozen. blade-structure stays separate: neither repo depends on the other.

The fork changed only what moving it needed: imports re-rooted from
`blade_structure_utilities.*_utilities` onto `blade_geometry_model`, and the
`tests` folders renamed. The files below were left out because they import
packages that stay in blade-structure (its shell mesh hand-off, its thin-wall model,
its section properties, its layup and loads, and `blade_structure_precompute`). They check the seam between this
geometry and the physics on top of it, so their counterparts belong in
`composite_laminate_model` and `anisotropic_section_model`, written against
those packages.

- `coordinate_frames/tests/test_cf002_sp001_exactness.py`
- `coordinate_frames/tests/test_frame_invariance_harness.py`
- `coordinate_frames/tests/test_ply_frame_transforms.py`
- `tests/section_geometry/test_blended_section_geometry.py`
- `section_geometry/examples/05_handoff_to_mesh.py`
- `section_geometry/tests/test_cap_feasibility_clamp.py`
- `section_geometry/tests/test_cap_ply_snap.py`
- `section_geometry/tests/test_cap_skin_bonds.py`
- `section_geometry/tests/test_cell_wall_loops.py`
- `section_geometry/tests/test_geometry_worker_facade.py`
- `section_geometry/tests/test_placement_anchoring.py`
- `section_geometry/tests/test_placement_wp6_metamorphic.py`
- `section_geometry/tests/test_shell_strip_topology.py`
- `section_geometry/tests/test_spar_cap_midline_offset.py`
- `section_geometry/visualisation/demos/cell_loops_thin_wall_demo.py`
- `section_geometry/visualisation/demos/downstream_consumers_demo.py`
- `section_geometry/visualisation/demos/geometry_to_mesh_flythrough.py`
- `section_geometry/visualisation/demos/quadtree_vs_shoelace_convergence_study.py`
- `section_geometry/visualisation/demos/section_properties_dashboard.py`
- `section_geometry/visualisation/demos/beyond_blades_demo.py`
- `section_geometry/visualisation/demos/design_vector_morph.py`
- `section_geometry/visualisation/demos/feature_flythroughs.py`
- `section_geometry/visualisation/demos/spanwise_web_anchors.py`
- `section_geometry/visualisation/demos/spanwise_web_anchors_shared.py`
- `section_geometry/visualisation/demos/spanwise_yz_plan_flyover.py`
