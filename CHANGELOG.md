# Changelog

## Unreleased

- `Stringer` on `BladeGeometry` is a free-end or through branch, optionally
  limited by `extent`. Inside that interval it appears on `ShellGraph.branches`;
  the ends are topology events and `forced_stations()`.

- `TransverseMember` on `BladeGeometry` places ribs, bulkheads and frames at a
  span station. `section_at` at that `z` attaches closed rib faces on
  `ShellGraph.faces` (cell-wall loops, or the skin loop when there are no
  cells). `forced_stations()` lists those `z` values; `topology_diff` records a
  rib event when the face is present on only one station. Goldens unchanged.

- Quadtree refinement is the array-based path only (`soa.py`). Object-model
  `refine.py` / `refine_multi.py` are gone; `build_area_mesh_2d` uses
  `build_shared_quadtree`.

- Live `docs/` no longer import or path into downstream packages; fork history
  stays in `docs/FORK.md`.

- Dual names warn until 0.6.0: `CellLoop` → `CellPolygon`, `station_at` →
  `interpolate`, `SectionGeometryRunner` → `SectionGeometry`,
  `web_x_positions_m` → `web_offsets`, `SpanWeb` → `ShearWeb`.

- `Airfoil` is the shape; `AirfoilSDF` is the implicit field via `Airfoil.to_sdf()`.
  NACA 4/5/6 vertex maths lives only in `sdf.shapes.naca`.

- Comments in `src/` and `tests/` use plain sentences instead of tracking codes
  (SG-/BL-/G-POLY/"Wave 0"); `docs/section/bug_log.md` still holds the history.

- CI: mypy on `src/implicit_geometry_kernel`, a consumer-name grep on `src/` and
  `tests/`, and `py.typed`. Branch protection on `main` requires the test jobs.

## 0.5.0

- Package version is `0.5.0`. `SpanWeb` remains an alias of `ShearWeb` until 0.6.0.

- Wall types register on `SECTION_BUILDERS` (`WallType` protocol). Label
  classification, midlines, junctions, cells and properties look up the
  registry instead of string prefixes. A dummy stiffener type registered only
  in tests flows through those stages. Goldens unchanged.

- Keyword aliases (`cap_height`, `web_thickness`, `twist_angle`) resolve on
  `SectionSpec.from_kwargs`; `kwarg_aliases.py` re-exports those helpers until
  0.6.0. Numbers unchanged.

- `section/api` folded into `section/spec` (moves only, geometry unchanged):
  runner, factory, layout resolution, results, kwarg aliases, and the worker
  live next to `SectionSpec`. `worker.py` is split into `build.py` (emit) and
  `prep.py` (cache key and sized assemble). `section.api` and its modules remain
  as aliases until 0.6.0.

- `SectionSpec`, `Gauges` and `SectionOptions` are the station geometry input
  (`section.spec`). `build_section_geometry(spec, options=...)` is the worker
  entry; the keyword form builds a spec and warns (removed in 0.6.0).
  `SectionGeometry` and `BladeGeometry.section_at` build a spec internally.
  Geometry unchanged.

- `blade/cut.py` split (moves only, geometry unchanged): wall/cell types stay in
  `cut.py`; the mid-surface graph, shell graph, partitions and topology diff live
  in `blade/shell_graph.py`. `blade.cut` still exports the graph names until 0.6.0.

- `section/placement` folded into `section/layout` (moves only, geometry unchanged):
  `anchoring.py`, `collision.py` and `handoff.py` now live beside the catalogue
  modules. `section.placement` and its modules remain as aliases until 0.6.0.

- `section/constructive` folded into its target homes (moves only, geometry unchanged):
  `export.py` is split into `section/boundary/walls.py` (skin, web and cap solids) and
  `section/boundary/cores.py` (cores, inserts, the combined set); `emitters.py` and
  `labels.py` move to `section/boundary/`; `properties.py` becomes
  `section/properties/polygons.py`. `section.constructive` and its modules remain as
  aliases until 0.6.0. `docs/REFACTOR_PLAN.md` status is updated: 5a/5b done, next is
  5c (`placement` into `layout`); Step↔Phase numbering is mapped in the status section.

- `section/` subpackages renamed to their target names (moves only, geometry
  unchanged): `section.props` is now `section.properties` and `section.recipes` is now
  `section.walls`. The old package and module paths remain as aliases until 0.6.0.

- `frames/` moved to its target tree (moves only, numbers unchanged): the 638-line
  `tensors.py` is now the package `frames/tensors/` with `inertia` (second moments and
  `structural_product_moment`), `moments`, `forces` and `thickness`.
  `frames.tensors` still exposes every former name.

- `integration/` moved to its target tree (moves only, geometry unchanged). The
  1,258-line `quadtree.py` and 566-line `quadrature.py` are now the package
  `integration/quadtree/`: `tree` (cells, meshes, roots), `refine` and `refine_multi`
  (object-model refinement), `soa` (array refinement, `build_shared_quadtree`),
  `marching` (zero contours, summaries), `moments` (exact cell moments) and `quadrature`
  (section properties). `integration.quadtree` still exposes every former name;
  `integration.quadrature` is an alias until 0.6.0.

- `sdf/` moved to its target tree (moves only, geometry unchanged): `sdf/nodes/`
  (`base`, `primitives`, `csg`, `bbox`, `tagged`, `span`, `lofted`), `sdf/kernels/`
  (`primitives`, `polygon_numba`, `grid`, `csg`, `transforms`) and `sdf/shapes/`
  (`airfoil`, `naca`). `sdf.ir` is split across `sdf.nodes`. The old module paths
  (`sdf.field`, `sdf.ir`, `sdf.span_nodes`, `sdf.lofted`, `sdf.primitives`, `sdf.grid`,
  `sdf.csg`, `sdf.transforms`, `sdf.airfoil`, `sdf.naca`) remain as aliases until 0.6.0;
  `sdf.nodes` is now the package and still exports the tagged nodes.
- `docs/REFACTOR_PLAN.md`: the restructure and exactness plan with a status table, so
  the work can be resumed from the repository alone.

- One Shapely gateway: `section.boundary` (module `section/boundary/polygons.py`) holds
  every polygon offset, Boolean clean-up, rotation and clip, and re-exports the Shapely
  types the section uses. No other section module imports Shapely (a test enforces it;
  `integration/shoelace.py` is the only other user). `section.constructive.helpers`,
  `section.constructive.offset_locus` and `section.midlines.offset_math` remain as
  aliases until 0.6.0. Geometry unchanged.
- One shape model: the CSG shapes in `sdf.ir` (`Circle`, `Polygon`, `Union`, `Offset`,
  `Rotate`, `Slab`, ...) are now `SDFNode` subclasses that evaluate themselves, carry
  feature ids through combinators and have a `bbox()` method. The compile step is gone:
  `compile_expr_to_node` returns its argument and `eval_via_nodes` equals
  `eval_expr_batch` (both kept until 0.6.0). `UnionNode`, `IntersectNode`,
  `SubtractNode`, `OffsetNode`, `RotateNode`, `ShellNode` and `SlabNode` are now the
  same classes as `Union`, `Intersect`, ... (positional arguments unchanged).
  `ExprSDF(expr, meta)` only stamps a feature id. Geometry unchanged.
- Sections are assembled through `structure`. `MultiCellSection` now plans its inputs
  (`section/assembly/plan.py`), runs member builders (`section/assembly/members.py`:
  skin, webs, caps, cores, inserts) and assembles the result; `MultiCellSection.structure`
  returns it and `registry=` adds builders. Geometry is unchanged (golden numbers
  identical). Caps are stacked wall layers under the skin.
- `structure`: `WallContribution(stacked=True)` for a layer beneath the wall; a partition
  now replaces the base wall inside its footprint (highest priority wins) instead of
  taking the larger thickness; every contribution takes an optional `source`.
- Safety net for the restructure: `tests/golden/` records section properties, midline
  strips and cell loops for every catalogue layout at 0, 15 and 35 degrees of twist, plus a
  lofted blade cut at five stations, and a snapshot of the public API. `tests/test_snapshots.py`
  fails if any of them move. Regenerate with `python scripts/snapshots.py --update` only for
  an intended change, and say why here.
- `scripts/check_consumer.py PATH` runs a downstream package's tests against this checkout.
- New `structure` package (geometry only): wall, branch, region and transverse
  contributions, a `MemberBuilder` extension point, `StructureAssembler`, and span laws
  (`ConstantLaw`, `TableLaw` with monotone PCHIP, and combinators). It assembles the wall
  thickness field, wall reference surface, cavity and a priority-owned material field,
  and reports the stations where members start or end. See `docs/structure/README.md`.
- New SDF nodes `ComplementNode`, `MaxNode`, `LawOffsetNode` and `SpanIntervalNode`;
  `FeatureKind` gains `core`, `stringer`, `transverse` and `bond`.
- Fixed: an opaque `CallableField` leaf now compares by its token and its callable. Comparing
  the token alone let the compile cache hand one section's clip function to a later section
  whose token matched (tokens built from `id()` repeat after garbage collection), so results
  could depend on what had been built earlier in the process.

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
