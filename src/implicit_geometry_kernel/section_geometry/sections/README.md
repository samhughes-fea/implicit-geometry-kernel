# sections — MultiCellSection assembly

↑ [`section_geometry`](../README.md)

Assembles a layout, aerofoil, and wall thicknesses into a `MultiCellSection`:
a labeled map of frozen wall recipes plus cached SDF evaluation. OML is the
aerofoil contour; IML is that contour offset by skin thickness. Twist is stored
for emit; solids are authored untwisted. Layout:
[`../layout/`](../layout/README.md). Recipes:
[`../recipes/`](../recipes/README.md). Thicknesses are geometric — materials are
assigned downstream (`section_bind`).

## Table of contents

1. [Organisation](#organisation)
2. [API](#api)
3. [Usage in the package](#usage-in-the-package)
4. [See also](#see-also)

---

## Organisation

| Path | Role |
|------|------|
| `multicell.py` | `MultiCellSection` host (all catalogue keys, N ≥ 0 webs) |
| `from_layout.py` | `multicell_from_layout` — `SectionLayout` + thicknesses → host |
| `labels.py` | `ComponentLabel` |
| `clip_algebra.py` | Public clip / half-plane helpers (wraps private `clips.py`) |
| `build_spar_caps.py` | `CapPlacementRule` → `CapWall` recipes → SDF IR |
| `web_layout.py` | Web anchor resolution |
| `cap_feasibility.py` | Skip caps that do not fit |

## API

| Name | File | Role |
|------|------|------|
| `MultiCellSection` | `multicell.py` | Assembled host (sole production host) |
| `multicell_from_layout` | `from_layout.py` | Layout → MultiCell (`SectionLayout` SoT) |

SDF IR classes (`OuterSkin`, `SparCap`, …) are **not** exported here. Import
from `recipes.sdf_ir` only inside projection / test code.

## Usage in the package

`SectionGeometryRunner` / `multicell_from_layout` are the production assemble
paths. Catalogue keys resolve to `SectionLayout` presets before assemble.
Caps are authored as `CapWall` recipes first; SDF components are a projection.

## See also

- [`../docs/reference/G_ROW_OWNERSHIP_MATRIX.md`](../docs/reference/G_ROW_OWNERSHIP_MATRIX.md)
- [`../layout/`](../layout/README.md)
- [`../recipes/`](../recipes/README.md)
- [`../sdf/`](../sdf/README.md)
