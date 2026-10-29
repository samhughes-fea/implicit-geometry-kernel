# types — public type aliases

↑ [`section_geometry`](../README.md)

Public type aliases for the package barrel. Definitions and behaviour live in
the owning modules; this home only re-exports stable names.

## Table of contents

1. [Organisation](#organisation)
2. [API](#api)
3. [See also](#see-also)

---

## Organisation

| File | Role |
|------|------|
| `__init__.py` | Public `__all__` of aliases |

## API

| Alias | Defined in | Role |
|-------|------------|------|
| `SectionGeometryResult` | [`../api/worker.py`](../api/worker.py) | Worker product |
| `SectionPropertiesMeta` | [`../api/worker.py`](../api/worker.py) | Properties-backend metadata |
| `SectionPropertiesReport` | [`../props_emit/export.py`](../props_emit/export.py) | Exportable report |
| `SDFCallable` | [`../layout/types.py`](../layout/types.py) | Protocol for `φ(x, y)` |

On `SectionGeometryResult`, `properties_map` is filled under `"properties"`;
`midline_strips` and `cell_wall_loops` under `"midlines"` / `"cell_wall_loops"`.
Midlines and cell-wall loops are already in the twisted frame. Full field tables:
[`../README.md` § Results](../README.md#results).

## See also

- [`../api/`](../api/README.md)
- Midline / `CellWallLoop` DTOs: [`../../section_contract/`](../../section_contract/README.md)
