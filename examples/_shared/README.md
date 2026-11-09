# _shared

Plot helpers for numbered examples (`draw_cut` solids, `draw_shell` conformal skin mid, wall colours, `save`).
`blade_viz` holds the blade figure tree (`cross_sections` / `csg_mesh` / `regions` / `region_cutaways` / `planform` / `full_blade`), suction-role filter, and planform sampling.
`debug_mesh.dump_debug_mesh` writes `mids.obj` / `mids.vtk` / `samples.npz` next to loft figures (optional `cutaway.*` stem).

↑ [examples](../README.md)

Not a runnable example. `scripts/run_examples.py` skips this folder.

## Figures

Numbered examples that call `save()` follow this contract (see section/02 and sdf/01):

- Claim in the title; what was checked in `note()`.
- Wall colours from `WALL_COLOURS` only (skin navy, cap red, web teal). No default colour cycle for walls.
- Units on axes. Equal aspect on every section-plane panel.
- If the caption says solids, the fill must be visible. Galleries may use *display gauges* (thicker than the printed manufacturing values) and must say so in `note()`.
- Integer feature ids use `meta_cmap` / `meta_legend`, never `tab10` on raw ids.
- Legend off the geometry (`legend_outside`) or a colour bar. No eight-entry legend on the trailing edge.
