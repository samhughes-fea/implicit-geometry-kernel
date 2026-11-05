# Changelog

## Unreleased

- Docs and the import-boundary test follow the split of composite-laminate-model into
  laminate-constitutive-model and shell-fem; the kernel imports neither.
- Seven new worked examples (02 to 08) with figures and READMEs: the layout catalogue,
  signed distances and CSG, section properties, the mid-surface and cells, web
  placement, airfoils and root transitions, and coordinate frames. Each checks the
  kernel against an independent reference. `examples/README.md` indexes them.
- Fixed: `frames.transform_second_moments_S_to_B` (and `_B_to_S`) gave the wrong sign to
  the Ixy term of `Ixx_B`, so Ixx + Iyy was not invariant under twist. It now matches
  integrating the rotated outline.
- Demo figures are written to `examples/demos/figures/`, committed and embedded in the
  demo READMEs, instead of a git-ignored `output/` folder or a folder inside the
  installed package. `save_figure` defaults to 150 dpi.

## 0.2.0 (unreleased)

The kernel is now the geometry of blade and wing cross-sections only, regrouped into
`blade`, `section`, `sdf`, `frames`, `integration`, `contract` and `plotting`. See
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
- Matplotlib is only needed for `plotting`; every other module imports without it.

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
