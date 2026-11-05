# shoelace_integration

`implicit_geometry_kernel.integration.shoelace` is a shared worker package in `blade_structure_utilities`.

## Purpose

- Provide reusable domain logic behind a stable package boundary.
- Keep cross-package imports routed through canonical worker paths.

## Public Interface

- Canonical imports use `implicit_geometry_kernel.integration.shoelace`.
- API barrel lives at `implicit_geometry_kernel.integration.shoelace.api`.
- Type/DTO barrel lives at `implicit_geometry_kernel.integration.shoelace.types`.

## Notes

- This README records the worker's canonical API/types surface and import-boundary contract.
- Expand this file with domain equations, ownership tables, and examples as work proceeds.
