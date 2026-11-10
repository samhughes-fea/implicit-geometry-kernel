# Quick start

Progressive call paths for one station. Overview: [`../README.md`](../README.md).
Recipes: [`caller_patterns.md`](caller_patterns.md).

`SectionGeometry.execute()` returns `SectionGeometryResult` or raises
`InfeasibleSectionGeometry`. Catalogue tokens are defined in the
[WX-Y-Z-NR taxonomy](../../blade/concepts/taxonomy.md).

↑ [`section`](../docs_index.md)

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
`../../../examples/section/teaching/01_first_station.py`.

```python
from implicit_geometry_kernel import Airfoil
from implicit_geometry_kernel.section import (
    SectionGeometry,
    SectionOptions,
    SectionSpec,
    resolve_system_type,
)
from implicit_geometry_kernel.plotting import plot_section

layout = resolve_system_type("2D-F-LC-0R")
airfoil = Airfoil.naca4("2412").to_sdf(chord=1.0)
spec = SectionSpec.from_kwargs(
    layout,
    airfoil,
    skin_thickness=0.004,
    spar_cap_thickness=0.012,
    shear_web_thickness=0.005,
    twist_rad=0.0,
)

result = SectionGeometry(
    spec,
    SectionOptions(outputs=("properties", "midlines", "cell_wall_loops")),
).execute().get_results()

plot_section(result.section)
```

`SectionGeometry(spec, options)` is the public call path. Keyword `layout=` /
`airfoil=` still work this minor and are removed in the next.
`build_section_geometry` is a deep import (`section.spec`) for tests. Gate
`"cell_wall_loops"` also fills midlines.

---

## 2. Read the result

Runnable twin (properties table):
`../../../examples/section/teaching/02_read_properties.py`.
Anatomy of labels on a figure:
[`anatomy_on_a_figure.md`](anatomy_on_a_figure.md).

| Field | When present | Use |
|-------|--------------|-----|
| `result.section` | always | Named planar solids (for plotting) |
| `result.properties_map` | `"properties"` gated | Per-region area / second moments |
| `result.midline_strips` | `"midlines"` or `"cell_wall_loops"` gated | Skin+web centrelines (conformal mesh) |
| `result.cell_wall_loops` | `"cell_wall_loops"` gated | Closed cell-wall loops |

```python
skin_area = result.properties_map["skin_suction_surface"]["area"]
for strip in result.midline_strips or ():
    print(strip.kind, strip.label, strip.midline_twisted_pitch_axis_blade_frame.shape)
for loop in result.cell_wall_loops or ():
    print(loop.id, loop.area_m2)
```

Field schemas: [Using the result](../docs_index.md#using-the-result) and
[Results](../docs_index.md#results) in the package README.

**Name map:** gate `"cell_wall_loops"` → field `cell_wall_loops`.

---

## 3. Pick a layout

Layout ids are strings in `SYSTEM_TYPE_KEYS` (120 keys: 40 cores × LC/LE/PA × `0R`).
A rib count is required; `"2D-F-LC-0R"` is none, `"2D-F-LC-10R"` is ten faces.
Section goldens iterate `SYSTEM_TYPE_LAYOUT_KEYS` (the 40 `*-LC-0R` cores). Examples:

| Need | Start with |
|------|------------|
| Twin-web torsion box | `"2D-F-LC-0R"` or `"2D-CN-LC-0R"` |
| Skin-only station | `"0A-LC-0R"` |
| Skin + single cap band | `"0B-LC-0R"` or `"1B-F-LC-0R"` |
| Three webs, continuous box | `"3D-F-LC-0R"` |

Browse every key as a figure:
[`../reference/structural_systems_gallery.md`](../reference/structural_systems_gallery.md).

Decode `[W][X]-[Y]-[Z]-[R]` (web count, cap family, orientation, LC/LE/PA, rib count):
[`taxonomy`](../../blade/concepts/taxonomy.md).

```python
from implicit_geometry_kernel.section import (
    SYSTEM_TYPE_KEYS,
    resolve_system_type,
)

assert "2D-F-LC-0R" in SYSTEM_TYPE_KEYS
spec = resolve_system_type("2D-F-LC-0R")
# spec.n_webs, spec.structural_family, spec.web_orientation, ...
```

Layouts `0A` / `0B` are valid zero-web hosts; shell-mesh handoff accepts
skin ± optional caps (no web junctions). See
[`ADR_ZERO_WEB_SHELL_MESH.md`](../reference/ADR_ZERO_WEB_SHELL_MESH.md).

---

## 4. Thickness-driven build

Thicknesses are metres on `SectionSpec`
(`skin_thickness`, `spar_cap_thickness`, `shear_web_thickness`).
Same layout, different stacks:

```python
spec = SectionSpec.from_kwargs(
    resolve_system_type("2D-F-LC-0R"),
    Airfoil.naca4("2412").to_sdf(chord=1.0),
    skin_thickness=0.006,
    spar_cap_thickness=0.020,
    shear_web_thickness=0.008,
)
result = SectionGeometry(
    spec,
    SectionOptions(outputs=("properties", "midlines")),
).execute().get_results()
```

Caller parameters and feasibility:
[`../spec/README.md`](../spec/README.md).

---

## 5. Sizing loop (resolve once)

Freeze which webs fit at a minimum skin, then vary thicknesses without
re-resolving so the design-variable set stays fixed. Full pattern:
[`caller_patterns.md`](caller_patterns.md#2-resolve-once-then-rebuild-sizing-loop).
Runnable twin:
`../../../examples/section/teaching/03_resolve_then_size.py`.
Feasibility notes: [`feasibility_and_drops.md`](feasibility_and_drops.md).

```python
frozen = SectionGeometry(
    SectionSpec.from_kwargs(
        layout, airfoil,
        skin_thickness=t_skin_min,
        spar_cap_thickness=0.0,
        shear_web_thickness=0.0,
    ),
).resolve().get_layout_resolution()

result = SectionGeometry(
    SectionSpec.from_kwargs(
        frozen.layout, airfoil,
        skin_thickness=skin_thickness,
        spar_cap_thickness=spar_cap_thickness,
        shear_web_thickness=shear_web_thickness,
    ),
    SectionOptions(outputs=("properties",)),
    resolution=frozen,
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
`../../../examples/section/teaching/04_spanwise_web_anchors.py`.

```python
result = SectionGeometry(
    SectionSpec.from_kwargs(
        layout, airfoil_at_station,
        skin_thickness=0.004,
        spar_cap_thickness=0.0,
        shear_web_thickness=0.0,
        web_offsets=(0.15, 0.50),  # metres from LE
    ),
    SectionOptions(outputs=("midlines",)),
).execute().get_results()
```

**Local fractions** (legacy): omit `web_offsets`. Positions become
`layout.web_chord_fracs × local_chord`.

---

## 7. What next

| Goal | Go to |
|------|--------|
| Runnable teaching `01`–`06` | [`../../../examples/section/teaching/`](../../../examples/section/teaching/) |
| More call recipes | [`caller_patterns.md`](caller_patterns.md) |
| Labels on a `"2D-F-LC-0R"` figure | [`anatomy_on_a_figure.md`](anatomy_on_a_figure.md) |
| Feasibility / dropped webs | [`feasibility_and_drops.md`](feasibility_and_drops.md) |
| Geometry → mesh walkthrough | [`geometry_to_mesh.md`](geometry_to_mesh.md) |
| All 40 layout figures | [`../reference/structural_systems_gallery.md`](../reference/structural_systems_gallery.md) |
| Layout id decode / cascade | [taxonomy](../../blade/concepts/taxonomy.md) |
| Module map | [`../README.md`](../docs_index.md#package-structure) |
| Run tests | `../tests/README.md` |
| Figure demos | `../../../examples/demos/section/` |

Handoff to meshing (requires midlines). Runnable twin:
[`../../../examples/section/teaching/05_midline_handoff.py`](../../../examples/section/teaching/05_midline_handoff.py).
Guide: [`geometry_to_mesh.md`](geometry_to_mesh.md).

```python
# Midlines are the mesh hand-off. A 3-D shell mesh is built downstream from
# result.midline_strips (and, on a loft, cut.shell_graph()).
for strip in result.midline_strips or ():
    _ = strip.label, strip.thickness_m, strip.midline_twisted_pitch_axis_blade_frame
```
