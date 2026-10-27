# tests

Unit and integration tests for `section_geometry`.

↑ [`section_geometry`](../README.md)

## Table of contents

1. [How to run](#how-to-run)
2. [What is covered (glance)](#what-is-covered-glance)

---

## How to run

From the repository root:

```bash
pytest blade_structure_utilities/geometry_utilities/section_geometry/tests -q
```

## What is covered (glance)

| Area | Examples |
|------|----------|
| Public facade / runner | `test_geometry_runner.py`, `test_geometry_worker_facade.py` |
| Layout resolve / cascade / web drop | `test_resolve_station_layout.py`, `test_layout_step_down_cascade.py`, `test_web_cavity_depth_drop.py`, `test_spanwise_web_anchoring.py` |
| Feasibility | `test_infeasible_section_geometry.py`, `test_cap_feasibility_clamp.py` |
| Assembly / families | `test_multicell.py`, `test_structural_family.py`, `test_skin_ss_ps_taxonomy.py` |
| SDF / CSG / IR | `test_airfoil.py`, `test_csg.py`, `test_csg_ir.py`, `test_primitives.py` |
| Constructive / props parity | `test_constructive_polygons.py`, `test_hybrid_parity.py`, `test_props_backend_parity.py` |
| Midlines / topology / cell-wall loops | `test_shell_strip_topology.py`, `test_cell_wall_loops.py`, `test_spar_cap_midline_offset.py` |
| Plots | `test_plot_section_b_frame.py`, `test_plot_from_mesh.py` |
