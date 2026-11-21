# 16 · Bottleneck splits

One parametric blade and one recorded quadtree station (NACA 2412, `2D-CN-LC-0R`, `backend="quadtree"`). The script
prints three splits. The figure is drawn and not saved.
NREL `h = 1.5` is a second block. It is skipped, with a one-line note, when
the reference-blades catalogue is absent. GEOS is required, as for the other
geometry examples.

```bash
python examples/blade/16_bottleneck_splits/timings.py
```

```text
shell_mesh 0.000s  figure 0.000s
shell_stations 0.000s  shell_strips 0.000s
stage_timings quadtree_phi=0.000s
sdf_polygon_count 0
```

The zeros above are the shape of the print. A run fills them. `shell_stations`
is the section-cut loop. `shell_strips` is the mapped-strip loop.
`quadtree_phi` is the φ batch inside quadtree refinement, separate from child
expansion. `sdf_polygon_count` is the profile counter for that station.
