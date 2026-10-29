# visualisation — plots and demos

↑ [`section_geometry`](../README.md)

Plot helpers and demos that sample and draw products from other modules.
Section figures often map twisted Y to the horizontal axis (LE on the right);
that remapping is display-only.

## Table of contents

1. [Organisation](#organisation)
2. [API](#api)
3. [See also](#see-also)

---

## Organisation

| Path | Role |
|------|------|
| `plot.py` | Public re-exports |
| `plot_section.py` | Contour / fill plots |
| `plot_sdf.py` | Raw SDF field |
| `plot_topology.py` | Midline overlays |
| `plot_subcomponents.py` | Per-label panels |
| `demos/` | Runnable scripts — [`demos/README.md`](demos/README.md) |
| `output/` | Generated figures |

## API

| Name | File | Role |
|------|------|------|
| `plot_section` | `plot_section.py` | Draw an assembled section |
| `plot_sdf_field` | `plot_sdf.py` | Draw a raw SDF field |

## See also

- Teaching scripts: [`../examples/`](../examples/README.md)
- [`demos/README.md`](demos/README.md)
- [`../../coordinate_frames/README.md`](../../coordinate_frames/README.md)
