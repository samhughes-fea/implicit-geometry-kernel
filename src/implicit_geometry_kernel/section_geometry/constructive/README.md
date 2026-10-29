# constructive — polygon path for area and inertia

↑ [`section_geometry`](../README.md)

For mature labels (skin, webs, caps, cores, LE/TE inserts), builds explicit
boundary polygons from the same recipes and integrates area / second moments
by shoelace. Residual φ remains available as a debug backend when an insert
recipe sets `emission_maturity='residual_phi'`. Parity:
[`../docs/reports/DUAL_REP_ACCEPTANCE.md`](../docs/reports/DUAL_REP_ACCEPTANCE.md). Packaging:
[`../props_emit/`](../props_emit/README.md).

Prefer `emit_wall_solid_polygons` / `is_wall_solid_label`. Wall solids are
geometric thickness regions (skin / web / spar cap), not composite layup —
material binding belongs to `section_bind` / `section_layup` / `section_schedule`.

## Table of contents

1. [Organisation](#organisation)
2. [API](#api)
3. [Usage in the package](#usage-in-the-package)
4. [See also](#see-also)

---

## Organisation

| File | Role |
|------|------|
| `export.py` | `emit_constructive_polygons` orchestration |
| `emitters.py` | Per-component polygon builders |
| `labels.py` | Constructive vs residual classification |
| `properties.py` | Shoelace moments from polygons |
| `helpers.py` | Shared geometry helpers |

## API

| Name | File | Role |
|------|------|------|
| `emit_constructive_polygons` | `export.py` | Mature labels → polygons |
| Label classifiers | `labels.py` | Constructive vs residual |
| `section_properties_from_polygons` | `properties.py` | Area and inertia from polygons |

## Usage in the package

```mermaid
flowchart LR
  MCS[MultiCellSection] --> Class[classify_labels]
  Class --> Poly[emit_polygons]
  Poly --> Shoe[shoelace]
  Shoe --> Props[props_emit]
```

## See also

- [`../props_emit/`](../props_emit/README.md)
- [`../docs/reports/DUAL_REP_ACCEPTANCE.md`](../docs/reports/DUAL_REP_ACCEPTANCE.md)
