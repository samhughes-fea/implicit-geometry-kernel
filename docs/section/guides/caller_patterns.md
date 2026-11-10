# Caller patterns

Station call recipes. Overview: [`../README.md`](../README.md). Progressive
guide: [`quickstart.md`](quickstart.md). Web anchoring quantities follow
[taxonomy Z](../../blade/concepts/taxonomy.md#z--web-anchoring).

↑ [`section_geometry`](../docs_index.md) · [`docs index`](../docs_index.md)

## Table of contents

1. [Simple build](#1-simple-build)
2. [Resolve once, then rebuild (sizing loop)](#2-resolve-once-then-rebuild-sizing-loop)
3. [Spanwise web positions vs local chord fractions](#3-spanwise-web-positions-vs-local-chord-fractions)

---

## 1. Simple build

One station, all optional outputs. Same as the package quick start.
Runnable twin: `../../../examples/section/teaching/01_first_station.py`.

```python
from implicit_geometry_kernel.blade import Airfoil
from implicit_geometry_kernel.section import (
    SectionGeometry,
    SectionOptions,
    SectionSpec,
    resolve_system_type,
)

spec = SectionSpec.from_kwargs(
    resolve_system_type("2D-F-LC-0R"),
    Airfoil.naca4("2412").to_sdf(chord=1.0),
    skin_thickness=0.004,
    spar_cap_thickness=0.012,
    shear_web_thickness=0.005,
)
result = SectionGeometry(
    spec,
    SectionOptions(outputs=("properties", "midlines", "cell_wall_loops")),
).execute().get_results()
```

## 2. Resolve once, then rebuild (sizing loop)

Freeze which webs fit at a minimum skin, then vary thicknesses without
re-resolving (so the design-variable set stays fixed). Runnable twin:
`../../../examples/section/teaching/03_resolve_then_size.py`.

```python
runner0 = SectionGeometry(
    SectionSpec.from_kwargs(
        layout, airfoil,
        skin_thickness=t_skin_min,
        spar_cap_thickness=0.0,
        shear_web_thickness=0.0,
    ),
)
frozen = runner0.resolve().get_layout_resolution()

for skin_thickness, spar_cap_thickness, shear_web_thickness in trials:
    result = SectionGeometry(
        SectionSpec.from_kwargs(
            frozen.layout, airfoil,
            skin_thickness=skin_thickness,
            spar_cap_thickness=spar_cap_thickness,
            shear_web_thickness=shear_web_thickness,
        ),
        SectionOptions(outputs=("properties",)),
        resolution=frozen,  # carries surviving web x/y
    ).execute().get_results()
```

If resolve cannot anchor any required web → `InfeasibleSectionGeometry`
(distinct from a thick trial).

## 3. Spanwise web positions vs local chord fractions

**Anchored (taper):** absolute offsets from LE, usually fixed from the root.
Webs that no longer fit are **dropped**, not slid. Runnable twin:
`../../../examples/section/teaching/04_spanwise_web_anchors.py`.
Guide: [`feasibility_and_drops.md`](feasibility_and_drops.md).

```python
result = SectionGeometry(
    SectionSpec.from_kwargs(
        layout, airfoil_at_station,
        skin_thickness=0.004,
        spar_cap_thickness=0.0,
        shear_web_thickness=0.0,
        web_offsets=(0.15, 0.50),  # metres from LE, fixed along span
    ),
    SectionOptions(outputs=("midlines",)),
).execute().get_results()
```

**Local fractions (legacy):** omit `web_offsets`. Positions become
`layout.web_chord_fracs × local_chord`. Every web scales with the section;
none is dropped for taper.

See also [`../spec/README.md`](../spec/README.md) (caller parameters) and
[`../layout/README.md`](../layout/README.md) (cascade).
