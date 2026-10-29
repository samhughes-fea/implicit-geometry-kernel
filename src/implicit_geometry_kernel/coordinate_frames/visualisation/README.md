# visualisation — plot demos and outputs

Operational suite for the coordinate-frames package: figure-producing demos and
generated PNGs. Frame mathematics live in
[`section_frames/`](../section_frames/README.md) and
[`shell_frames/`](../shell_frames/README.md) (see [`../README.md`](../README.md)).

## Purpose / ownership

| Owns | Does **not** own |
|------|------------------|
| Runnable demos under `demos/` and PNGs under `output/` | S↔B / E↔M transform formulas |
| Demo catalog and run commands | Product FE fringe / MITC4 dashboards |

Demos **must** import production transforms from the package barrel — never
reimplement NACA↔S / S↔B / plot-axis maths locally.

## Directory map

| Path | Role |
|------|------|
| `demos/` | Runnable figure demos |
| `output/` | Generated figures |
| `__init__.py` | Convenience re-export of section plot helpers |

## Demos

| Script | What it proves | Output |
|--------|----------------|--------|
| [`demos/plot_structural_frames.py`](demos/plot_structural_frames.py) | Shared B/S/S′ storage + plot contract (LE on RHS, \(+Y\)→TE); panel order B / S / S′ | `output/structural_frames_B_S_Sprime.png` |

```bash
python -m blade_structure_utilities.geometry_utilities.coordinate_frames.visualisation.demos.plot_structural_frames
```

## See also

- Section frames: [`../section_frames/README.md`](../section_frames/README.md)
- Shell frames: [`../shell_frames/README.md`](../shell_frames/README.md)
- Plot contract API: `section_frames.plot` (also on the package barrel)
