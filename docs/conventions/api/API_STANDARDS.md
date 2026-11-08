# API standards

**Audience: internal.** Worker-runner notes for maintainers.

↑ [conventions](../README.md)

Public station geometry is `SectionGeometry` on `implicit_geometry_kernel.section`.
`build_section_geometry` is a deep import (`section.spec`) for tests.

## Worker runners

A runner owns one execute path:

1. Construct with a spec (layout, airfoil, gauges, position, `outputs` gates).
2. `execute()` builds solids and gated products.
3. `get_results()` returns `SectionGeometryResult`, or the call raises
   `InfeasibleSectionGeometry`.

`resolve()` is the layout-only pass (web drops) without assembling solids.

Gates are `"properties"`, `"midlines"`, `"cell_wall_loops"` — products, not use-case
modes. See [docs/section/spec](../../section/spec/README.md).

Root barrel exports loft + catalogue keys. Import `SectionGeometry` from
`.section`, contract types from `.contract`.
