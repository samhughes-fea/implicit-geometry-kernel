# tests/section

Unit and integration tests for `implicit_geometry_kernel.section`.

↑ [tests](../README.md) · [docs](../../docs/section/README.md) ·
[examples](../../examples/section/README.md)

```bash
pytest tests/section -q
```

| Area | Examples |
|------|----------|
| Public facade / runner | `test_geometry_runner.py` |
| Layout resolve / cascade / web drop | `test_resolve_station_layout.py`, `test_spanwise_web_anchoring.py` |
| Feasibility | `test_infeasible_section_geometry.py` |
| Assembly / families | `test_multicell.py`, `test_structural_family.py` |
| Constructive / props parity | `test_constructive_polygons.py`, `test_hybrid_parity.py` |
| Midlines / cell-wall loops | `test_shell_strip_topology.py`, `test_cell_wall_loops.py` |
| Plots | `test_plot_section_b_frame.py` |
