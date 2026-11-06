# Where this geometry came from

`section_geometry`, `section_contract`, `coordinate_frames`, `area_2d` and
`shoelace_integration` are a one-time fork, with history, of the same packages in
[samhughes-fea/blade-structure](https://github.com/samhughes-fea/blade-structure)
at commit `51ba12f0` (27 October 2026), where `section_geometry` is complete and
frozen. blade-structure stays separate: neither repo depends on the other.

The fork changed only what moving it needed: imports re-rooted from
`blade_structure_utilities.*_utilities` onto `implicit_geometry_kernel`, and the
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

## After the split (29 October 2026)

The package moved from `blade_geometry_model` in composite-laminate-fea to its own
repository as `implicit_geometry_kernel`; only import lines changed. The fork's
package-level test folders (`tests/coordinate_frames`, `tests/section_geometry`) moved
to the repository's top-level `tests/`.

## Reorganisation (0.2.0, November 2026)

The kernel's role was fixed as the geometry of blade and wing cross-sections, with no
materials, plies, stiffness, loads or solvers, and the package was regrouped around
that. The file lists above use the names from before this change.

| Before | After |
|---|---|
| `blade_geometry.py`, `loft.py` | `blade/` (`airfoil.py`, `stations.py`, `loft.py`, `cut.py`) |
| `section_geometry/` | `section/` |
| `section_geometry/sdf/` | `sdf/` |
| `section_geometry/sdf/cache.py` | `section/eval_cache.py` |
| `section_geometry/profiling.py` | `profiling.py` |
| `section_geometry/visualisation/`, `coordinate_frames/plot.py` | `plotting/` (sibling of the kernel, not under `src/`) |
| `coordinate_frames/` | `frames/` |
| `section_contract/` | `contract/` |
| `area_2d/`, `shoelace_integration/` | `integration/` (shoelace in `integration/shoelace.py`) |
| each subpackage's `tests/`, `examples/`, `docs/` | top-level `tests/`, `examples/`, `docs/` |

Deleted: the `api/` and `types/` re-export barrels, `section_geometry/topology/` and other
compatibility aliases, the placeholder `example_plot.py` demos, functions with no
caller in any of the four repositories, and the code that belongs to the physics on top
of the geometry: ply stress and strain rotations (`shell_frames/ply.py`, duplicated in
composite-laminate-model), the shell-solver load channels (`shell_frames/load_channels.py`),
the beam stiffness transform (`section_frames/stiffness.py`) and `BeamResultantVector`
with its load-order constants.

## 0.3.0

0.3.0 removed the deprecated import paths (`implicit_geometry_kernel.blade_geometry`,
`implicit_geometry_kernel.section_geometry`, and the section shim packages). Import
`blade`, `section`, `sdf`, `frames`, `integration`, and `contract` from
`implicit_geometry_kernel` directly.
