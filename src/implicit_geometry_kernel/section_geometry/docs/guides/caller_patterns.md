# Caller patterns

Call recipes. Overview: [`../README.md`](../../README.md). Progressive guide:
[`quickstart.md`](quickstart.md).

↑ [`section_geometry`](../../README.md) · [`docs index`](../README.md)

## Table of contents

1. [Simple build](#1-simple-build)
2. [Resolve once, then rebuild (sizing loop)](#2-resolve-once-then-rebuild-sizing-loop)
3. [Spanwise web positions vs local chord fractions](#3-spanwise-web-positions-vs-local-chord-fractions)

---

## 1. Simple build

One station, all optional outputs. Same as the package quick start.
Runnable twin: [`../examples/01_first_station.py`](../../examples/01_first_station.py).

```python
from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    SectionGeometryRunner,
    resolve_system_type,
)

result = SectionGeometryRunner(
    layout=resolve_system_type("2D-F"),
    airfoil=AirfoilSDF.from_naca("2412", chord=1.0),
    skin_thickness=0.004,
    spar_cap_thickness=0.012,
    shear_web_thickness=0.005,
    outputs=("properties", "midlines", "cell_wall_loops"),
).execute().get_results()
```

## 2. Resolve once, then rebuild (sizing loop)

Freeze which webs fit at a minimum skin, then vary thicknesses without
re-resolving (so the design-variable set stays fixed). Runnable twin:
[`../examples/03_resolve_then_size.py`](../../examples/03_resolve_then_size.py).

```python
runner0 = SectionGeometryRunner(
    layout=layout,
    airfoil=airfoil,
    skin_thickness=t_skin_min,
    spar_cap_thickness=0.0,
    shear_web_thickness=0.0,
)
frozen = runner0.resolve().get_layout_resolution()

for skin_thickness, spar_cap_thickness, shear_web_thickness in trials:
    result = SectionGeometryRunner(
        layout=frozen.layout,
        airfoil=airfoil,
        skin_thickness=skin_thickness,
        spar_cap_thickness=spar_cap_thickness,
        shear_web_thickness=shear_web_thickness,
        resolution=frozen,  # carries surviving web x/y
        outputs=("properties",),
    ).execute().get_results()
```

If resolve cannot anchor any required web → `InfeasibleSectionGeometry`
(distinct from a thick trial).

## 3. Spanwise web positions vs local chord fractions

**Anchored (taper):** absolute offsets from LE, usually fixed from the root.
Webs that no longer fit are **dropped**, not slid. Runnable twin:
[`../examples/04_spanwise_web_anchors.py`](../../examples/04_spanwise_web_anchors.py).
Guide: [`feasibility_and_drops.md`](feasibility_and_drops.md).

```python
result = SectionGeometryRunner(
    layout=layout,
    airfoil=airfoil_at_station,
    skin_thickness=0.004,
    web_x_positions_m=(0.15, 0.50),  # metres from LE, fixed along span
    outputs=("midlines",),
).execute().get_results()
```

**Local fractions (legacy):** omit `web_x_positions_m`. Positions become
`layout.web_chord_fracs × local_chord`. Every web scales with the section;
none is dropped for taper.

See also [`../api/README.md`](../../api/README.md) (caller parameters) and
[`../layout/README.md`](../../layout/README.md) (cascade).
