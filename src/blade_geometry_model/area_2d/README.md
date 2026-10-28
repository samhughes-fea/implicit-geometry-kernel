# area_2d

`blade_structure_utilities.mesh_utilities.area_2d` is a shared worker package in `blade_structure_utilities`.

## Purpose

- Provide reusable domain logic behind a stable package boundary.
- Keep cross-package imports routed through canonical worker paths.

## Public Interface

- Canonical imports use `blade_structure_utilities.mesh_utilities.area_2d`.
- API barrel lives at `blade_structure_utilities.mesh_utilities.area_2d.api`.
- Type/DTO barrel lives at `blade_structure_utilities.mesh_utilities.area_2d.types`.

## Notes

- This README records the worker's canonical API/types surface and import-boundary contract.
- Expand this file with domain equations, ownership tables, and examples as work proceeds.
