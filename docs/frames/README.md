# Coordinate frames

↑ [docs](../README.md) · [src](../../src/implicit_geometry_kernel/frames/README.md) ·
[examples](../../examples/frames/README.md) · [tests](../../tests/frames/README.md)

The `frames` package is the single home of the section frames R, B, S, S′ and C,
their short tags (`FRAME_*`) and every transform between them. The definitions and
derivations are in **[R1 · Coordinate frames and twist](../r/R1_frames.md)**; this
page only says which array is in which frame. Do not reimplement the twist rotation,
the pitch-axis shift or the moment polarity outside this package.

**Import:** `from implicit_geometry_kernel.frames import …`; plot helpers
`section_plot_axes` and `setup_section_axes` are in
`implicit_geometry_kernel.plotting.frames`.

## Which array

| Object | Frame |
|---|---|
| `Station` scalars and aerofoil vertices | Aerofoil coordinates $(x_c, y_t)$ |
| `SectionGeometry` solids, `properties_map` | Authoring frame near S (LE-origin rotation when twisted) |
| `SectionCut.walls`, `CellPolygon.points` | B in swapped $(x, y)$ order, twist $\tau=-\theta$, (R1.11) |
| `ShellMidlineStrip`, `CellWallLoop` | B (`*_twisted_pitch_axis_blade_frame`) |
| `shell_graph` runs | B |

Known and fixed defects: [bug log](bug_log.md). Detail of the section-level
frame handling: [section_frames](section_frames/README.md); programme inventory:
[FRAMES_INVENTORY](FRAMES_INVENTORY.md).
