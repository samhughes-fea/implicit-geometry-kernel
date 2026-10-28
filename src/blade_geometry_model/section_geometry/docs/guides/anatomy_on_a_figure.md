# Anatomy on a figure

Map a twin-web `"2D-F"` plot to `SectionGeometryResult` labels. Overview:
[`../README.md`](../../README.md). Runnable twins:
[`../examples/01_first_station.py`](../../examples/01_first_station.py),
[`../examples/02_read_properties.py`](../../examples/02_read_properties.py).

↑ [`section_geometry`](../../README.md) · [`docs index`](../README.md)

## Table of contents

1. [Figure](#figure)
2. [Wall solids (`properties_map` / section labels)](#wall-solids)
3. [Midline strips](#midline-strips)
4. [Cell-wall loops](#cell-wall-loops)
5. [See also](#see-also)

---

## Figure

Layout `"2D-F"`: two shear webs, family-D continuous box caps, flapwise webs.
Demo thicknesses (NACA 2412, chord 1 m). Legend names the wall solids listed
below.

![2D-F labelled](../../visualisation/output/anatomy_2d_f/2D-F_labelled.png)

Generated, not checked in. Regenerate with:

```bash
python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/anatomy_2d_f_demo.py
```

Legend-less gallery tiles (all 40 keys) still come from
[`structural_family_grids_demo.py`](../../visualisation/demos/structural_family_grids_demo.py);
see [`ARTIFACTS.md`](../../visualisation/demos/ARTIFACTS.md).

Chord runs roughly horizontal in the section plot; LE is toward the right in
the structural plot contract (see
[`coordinate_frames`](../../../coordinate_frames/README.md)).

---

<a id="wall-solids"></a>
## Wall solids (`properties_map` / section labels)

Thicknesses are geometric wall gauges only (material-bind agnostic). Typical
`"2D-F"` keys:

| Label | What you see on the figure |
|-------|----------------------------|
| `skin_suction_surface` | Upper outer skin wall |
| `skin_pressure_surface` | Lower outer skin wall |
| `shear_web_0`, `shear_web_1` | LE-ward then TE-ward web solids |
| `spar_cap_suction_surface` | Upper box-cap band between the webs |
| `spar_cap_pressure_surface` | Lower box-cap band |
| `sandwich_core_*` | Optional core bays (when present in the build) |

```python
skin_area = result.properties_map["skin_suction_surface"]["area"]
```

Schemas: package README [Results](../../README.md#results).

---

## Midline strips

Gate `"midlines"` (or `"cell_wall_loops"`). Each `ShellMidlineStrip` is a
centreline polyline for meshing:

| `kind` | Typical `label` pattern |
|--------|-------------------------|
| `"skin"` | `skin_suction_surface:leading_edge`, `:bay_0`, `:trailing_edge` (and pressure-side twins) |
| `"web"` | `shear_web_0`, `shear_web_1` |
| `"cap"` | spar-cap suction / pressure labels |

Coordinates live in the twisted blade frame
(`midline_twisted_pitch_axis_blade_frame`).

---

## Cell-wall loops

Gate `"cell_wall_loops"`. Closed loops of the structural cavities, LE → TE:

| `id` | Role on `"2D-F"` |
|------|------------------|
| `0` | Leading-edge cell (ahead of the first web) |
| `1` | Box cell between the two webs |
| `2` | Trailing-edge cell (aft of the second web) |

`area_m2` is the shoelace area of the twisted loop, not a wall solid area.

---

## See also

- Gallery of all keys: [`../reference/structural_systems_gallery.md`](../reference/structural_systems_gallery.md)
- Taxonomy decode: [`../reference/system_type_xyz_taxonomy.md`](../reference/system_type_xyz_taxonomy.md)
- Midlines module: [`../midlines/`](../../midlines/README.md)
- Cell-wall loops module: [`../cell_wall_loops/`](../../cell_wall_loops/README.md)
