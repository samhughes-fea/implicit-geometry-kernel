# Quick start

Progressive call paths. Overview: [`../README.md`](../../README.md). Recipes:
[`caller_patterns.md`](caller_patterns.md).

↑ [`section_geometry`](../../README.md)

Returns `SectionGeometryResult` or raises (e.g. `InfeasibleSectionGeometry`).

## Table of contents

1. [First call](#1-first-call)
2. [Read the result](#2-read-the-result)
3. [Pick a layout](#3-pick-a-layout)
4. [Thickness-driven build](#4-thickness-driven-build)
5. [Sizing loop (resolve once)](#5-sizing-loop-resolve-once)
6. [Spanwise web positions](#6-spanwise-web-positions)
7. [What next](#7-what-next)

---

## 1. First call

One station. Runnable twin:
[`../examples/01_first_station.py`](../../examples/01_first_station.py).

```python
from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    SectionGeometryRunner,
    plot_section,
    resolve_system_type,
)

layout = resolve_system_type("2D-F")
airfoil = AirfoilSDF.from_naca("2412", chord=1.0)

result = SectionGeometryRunner(
    layout=layout,
    airfoil=airfoil,
    skin_thickness=0.004,
    spar_cap_thickness=0.012,
    shear_web_thickness=0.005,
    twist_rad=0.0,
    outputs=("properties", "midlines", "cell_wall_loops"),
).execute().get_results()

plot_section(result.section)
```

`SectionGeometryRunner` is the public call path; `build_section_geometry` is
the test surface. Gate `"cell_wall_loops"` also fills midlines.

---

## 2. Read the result

Runnable twin (properties table):
[`../examples/02_read_properties.py`](../../examples/02_read_properties.py).
Anatomy of labels on a figure:
[`anatomy_on_a_figure.md`](anatomy_on_a_figure.md).

| Field | When present | Use |
|-------|--------------|-----|
| `result.section` | always | Named planar solids (for plotting) |
| `result.properties_map` | `"properties"` gated | Per-region area / second moments |
| `result.midline_strips` | `"midlines"` or `"cell_wall_loops"` gated | Centreline polylines for meshing |
| `result.cell_wall_loops` | `"cell_wall_loops"` gated | Closed cell-wall loops |

```python
skin_area = result.properties_map["skin_suction_surface"]["area"]
for strip in result.midline_strips or ():
    print(strip.kind, strip.label, strip.midline_twisted_pitch_axis_blade_frame.shape)
for loop in result.cell_wall_loops or ():
    print(loop.id, loop.area_m2)
```

Field schemas: [Using the result](../../README.md#using-the-result) and
[Results](../../README.md#results) in the package README.

**Name map:** gate `"cell_wall_loops"` → field `cell_wall_loops`.

---

## 3. Pick a layout

Layout ids are strings in `SYSTEM_TYPE_KEYS` (40 keys). Examples:

| Need | Start with |
|------|------------|
| Twin-web torsion box | `"2D-F"` or `"2D-CN"` |
| Skin-only station | `"0A"` |
| Skin + single cap band | `"0B"` or `"1B-F"` |
| Three webs, continuous box | `"3D-F"` |

Browse every key as a figure:
[`../reference/structural_systems_gallery.md`](../reference/structural_systems_gallery.md).

Decode `{X}{Y}-{Z}` (web count, cap family, orientation):
[`../reference/system_type_xyz_taxonomy.md`](../reference/system_type_xyz_taxonomy.md).

```python
from blade_structure_utilities.geometry_utilities.section_geometry import (
    SYSTEM_TYPE_KEYS,
    resolve_system_type,
)

assert "2D-F" in SYSTEM_TYPE_KEYS
spec = resolve_system_type("2D-F")
# spec.n_webs, spec.structural_family, spec.web_orientation, ...
```

Layouts `0A` / `0B` are valid zero-web hosts; shell-mesh handoff accepts
skin ± optional caps (no web junctions). See
[`ADR_ZERO_WEB_SHELL_MESH.md`](../reference/ADR_ZERO_WEB_SHELL_MESH.md).

---

## 4. Thickness-driven build

Thicknesses are metres on the runner
(`skin_thickness`, `spar_cap_thickness`, `shear_web_thickness`).
Same layout, different stacks:

```python
result = SectionGeometryRunner(
    layout=resolve_system_type("2D-F"),
    airfoil=AirfoilSDF.from_naca("2412", chord=1.0),
    skin_thickness=0.006,
    spar_cap_thickness=0.020,
    shear_web_thickness=0.008,
    outputs=("properties", "midlines"),
).execute().get_results()
```

Caller parameters and feasibility:
[`../api/README.md`](../../api/README.md).

---

## 5. Sizing loop (resolve once)

Freeze which webs fit at a minimum skin, then vary thicknesses without
re-resolving so the design-variable set stays fixed. Full pattern:
[`caller_patterns.md`](caller_patterns.md#2-resolve-once-then-rebuild-sizing-loop).
Runnable twin:
[`../examples/03_resolve_then_size.py`](../../examples/03_resolve_then_size.py).
Feasibility notes: [`feasibility_and_drops.md`](feasibility_and_drops.md).

```python
frozen = SectionGeometryRunner(
    layout=layout,
    airfoil=airfoil,
    skin_thickness=t_skin_min,
    spar_cap_thickness=0.0,
    shear_web_thickness=0.0,
).resolve().get_layout_resolution()

result = SectionGeometryRunner(
    layout=frozen.layout,
    airfoil=airfoil,
    skin_thickness=skin_thickness,
    spar_cap_thickness=spar_cap_thickness,
    shear_web_thickness=shear_web_thickness,
    resolution=frozen,
    outputs=("properties",),
).execute().get_results()
```

If resolve cannot anchor any required web, it raises
`InfeasibleSectionGeometry`.

---

## 6. Spanwise web positions

**Anchored** (preferred for tapering blades): pass absolute offsets from the
leading edge. Webs that no longer fit are **dropped**, not slid. Details:
[`caller_patterns.md`](caller_patterns.md#3-spanwise-web-positions-vs-local-chord-fractions).
Runnable twin:
[`../examples/04_spanwise_web_anchors.py`](../../examples/04_spanwise_web_anchors.py).

```python
result = SectionGeometryRunner(
    layout=layout,
    airfoil=airfoil_at_station,
    skin_thickness=0.004,
    web_x_positions_m=(0.15, 0.50),  # metres from LE
    outputs=("midlines",),
).execute().get_results()
```

**Local fractions** (legacy): omit `web_x_positions_m`. Positions become
`layout.web_chord_fracs × local_chord`.

---

## 7. What next

| Goal | Go to |
|------|--------|
| Runnable examples `01`–`05` | [`../examples/`](../../examples/README.md) |
| More call recipes | [`caller_patterns.md`](caller_patterns.md) |
| Labels on a `"2D-F"` figure | [`anatomy_on_a_figure.md`](anatomy_on_a_figure.md) |
| Feasibility / dropped webs | [`feasibility_and_drops.md`](feasibility_and_drops.md) |
| Geometry → mesh walkthrough | [`geometry_to_mesh.md`](geometry_to_mesh.md) |
| All 40 layout figures | [`../reference/structural_systems_gallery.md`](../reference/structural_systems_gallery.md) |
| Layout id decode / cascade | [`../reference/system_type_xyz_taxonomy.md`](../reference/system_type_xyz_taxonomy.md) |
| Module map | [`../README.md`](../../README.md#package-structure) |
| Run tests | [`../tests/README.md`](../../tests/README.md) |
| Figure demos | [`../visualisation/demos/`](../../visualisation/demos/README.md) |

Handoff to meshing (requires midlines). Runnable twin:
[`../examples/05_handoff_to_mesh.py`](../../examples/05_handoff_to_mesh.py).
Guide: [`geometry_to_mesh.md`](geometry_to_mesh.md).

```python
from blade_structure_utilities.mesh_utilities.section_mesh import (
    shell_mesh_inputs_from_geometry_result,
)
shell_inputs = shell_mesh_inputs_from_geometry_result(result)
```
