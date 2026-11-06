# Changelog

## Unreleased

- Fixed: `Airfoil.naca4(code, n)` again samples `n` points per surface (0.4.0 passed `n`
  as the total, halving the density and opening a small gap at the leading edge).
- `ShearWeb` is the public name of the blade-station web again (top level and `blade`), so
  shell-fem and anisotropic-section-model import unchanged. `SpanWeb` (the 0.4.0 name) is an
  alias, removed in 0.5.0. The section-level SDF class in `section.recipes.shear_web` that
  shared the name is now `ShearWebComponent`.
- Lint: two unused imports in tests (CI was red on 0.4.0).

## 0.4.0

- README names the public contract: two input tracks (parametric vs windIO) into
  one assembler; `BladeGeometry` interpolates along span then calls
  `SectionGeometry` for one station; three products (`properties_map`,
  `ShellMidlineStrip` mesh hand-off, `CellWallLoop`).
- `SectionGeometryRunner` is an alias of `SectionGeometry`.
  `BladeGeometry.interpolate` is the scalar/foil interpolator (`station_at`
  remains an alias). Loft cell interiors are `CellPolygon` (`CellLoop` alias).
  Station metres-from-LE webs accept `web_offsets` (same quantity as
  `StationLayoutResolution.web_x_positions_m`).
- `Airfoil.naca4` and `AirfoilSDF.from_naca` share `naca_four_digit_vertices`.
  Prefer `Airfoil.to_sdf()`. `sdf` still does not import `blade`.
- `extract/` is documented (`docs/extract/`); off the kernel barrel. Node zoo
  imports from `sdf.nodes`. G-row documents `placement/` vs `CapPlacementRule`.
- Matplotlib helpers moved out of the kernel to the sibling package `plotting/`
  (`import plotting`, extra `[plots]`).
- Docs and the import-boundary test follow the split of composite-laminate-model into
  laminate-constitutive-model and shell-fem; the kernel imports neither.
- Worked examples 02–09 with figures and READMEs. `examples/README.md` indexes them.
- Fixed: `frames.transform_second_moments_S_to_B` (and `_B_to_S`) gave the wrong sign to
  the Ixy term of `Ixx_B`, so Ixx + Iyy was not invariant under twist.

## 0.3.0

Removed the deprecated import paths listed in 0.2.x (`blade_geometry`,
`section_geometry`, `frames.section_frames`, and the section shims
`sections` / `junction_preparation` / `cell_wall_loops` / `props_emit` /
`recipes.sdf_ir`). Offset math lives under `section.midlines`; G-POLY may
import that leaf. Public docs prefer `Airfoil.to_sdf()` and `SpanWeb`.

## 0.2.0

The kernel is now the geometry of blade and wing cross-sections only, regrouped into
`blade`, `section`, `sdf`, `frames`, `integration` and `contract`. See
[docs/FORK.md](docs/FORK.md) for the full map from the old names.

- `SYSTEM_TYPE_KEYS`, `resolve_system_type` and `MidSurface` are exported from the top
  level.
- Breaking: subpackage import paths change (`section_geometry` -> `section`,
  `coordinate_frames` -> `frames`, `section_contract` -> `contract`, `area_2d` and
  `shoelace_integration` -> `integration`, `section_geometry.sdf` -> `sdf`,
  `section_geometry.visualisation` -> `plotting`). `blade_geometry` and
  `section_geometry.layout.system_layout` keep working with a `DeprecationWarning`
  until 0.3.0.
- Removed: ply stress and strain rotations, shell-solver load channels, the beam
  stiffness transform, `BeamResultantVector` and its load-order constants (they belong
  to composite-laminate-model and anisotropic-section-model), the `api/` and `types/`
  barrels, the `topology` aliases, and code with no caller.
- Tests, examples, demos and docs moved out of the installed package; the wheel goes
  from 331 files to about 120.
- Matplotlib is only needed for the sibling `plotting` package; the kernel imports without it.

## 0.1.0 (29 October 2026)

- Split out of [composite-laminate-fea](https://github.com/samhughes-fea/composite-laminate-fea)
  with its history. The package was `blade_geometry_model` there and is now
  `implicit_geometry_kernel`.
- `BladeGeometry` lofts any of the 40 layouts along the span; `section_at(z)` returns
  the wall midlines, thicknesses, labels and cells; `SectionCut.mid_surface()` gives the
  wall graph that the laminate and section models build on.
- Flapwise (`-F`) webs stay normal to the rotor plane under twist.
- Forked from blade-structure: `section_geometry`, `section_contract`,
  `coordinate_frames`, `area_2d` and `shoelace_integration`.
