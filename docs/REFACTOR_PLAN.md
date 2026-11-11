# Refactor and exactness plan

**Audience: internal.** Maintainer notes; not part of the library manual.

This file is the working plan for restructuring the kernel and closing its exactness
gaps. It is written so that anyone (a person or an automated workflow) can pick the work
up from the **Status** section without other context. Update the status table in the
same PR that finishes a step.

## Status (resume here)

The architecture programme is **closed**. Steps 1–8 are done; step 9 was declined.
Standing work is product/API/docs cleanup, demo ports onto the public entry points,
and parked consumer items — not a new numbered step.

### 0.9.x plan (after 0.9.0)

The next programme: consumers pin the kernel and the kernel ships the analysis
shell mesh. Releases stay on 0.9.x (0.9.1, 0.9.2, ...); the maintainer will not
cut 1.0.0 while the public API is still being worked out. Kernel PRs:

| Id | What | State |
|----|------|-------|
| K1 | 0.9.0 release notes, `docs/MIGRATING.md`, 0.2 compat check | Done (#52, 0.9.1) |
| K2 | Reference-catalogue tests (`-m reference`), stable cap names at family-C stations, `pyyaml` in `[dev]` | Done (#52, 0.9.1) |
| K3 | No consumer names in the public repo; report renamed `SHELL_FE_REPLACEMENT_GAPS.md` | Done (#52, 0.9.1) |
| K4 | Pinned `ruff`/`mypy`, weekly `tool-drift` workflow | Done (#52, 0.9.1) |
| K6 | Tip closure on every blade: closes the skin shell at z/L = 1 with the polygon mesher, skin partition; `nR` counts interior ribs only | In review (0.9.2) |
| K7 | Mesh product: skin/web quads, polygonal ribs and tip, arrays with partition and thickness, no Gmsh | After K6 |
| K8 | Leading-edge and spanwise refinement controls | After K7 |
| K5 | Trimmed top-level API (deprecate first), internal docs to `docs/internal/`; still 0.9.x | After consumers move |

Maintainer decisions (10 Nov): the kernel produces the analysis shell mesh and
consumers use it without remeshing (Gmsh only for comparison in consumers);
skin and web mid-surfaces are quad patches, ribs and the tip closure polygonal
meshes; the tip closure is part of the skin shell, not a rib.

Status keys are **Steps 1–9**. Phase numbers in the body are historical names for
the same work; they are not extra steps.

| Step | Phase names in the body | What |
|------|-------------------------|------|
| 1 | — | Stabilise |
| 2 | — | Safety net |
| 3 | — | `structure/` member model |
| 4 | — | Sections assembled through `structure/` |
| 5 | Phase 1 (done), Phase 2 / 2d (after the remaining moves), Phase 3d | Target tree, then `SectionSpec` / `SectionOptions` |
| 6 | — | Seams with materials and meshing |
| 7 | Phase 5 | Completeness |
| 8 | Phase 4 | Exactness |
| 9 | — | Optional separate blade-model package |
| — | Phase 3a–c, 3e | Naming and plots; 3e is plots into `[plots]` |

| Step | What | State |
|------|------|-------|
| 1 | Stabilise: lint fix, naming fix, data clean-up | Done (PR #4, #5). Branch protection, consumer-name grep and mypy on `src/` added. `v0.4.1` was skipped; tagged `v0.5.0` |
| 2 | Safety net: golden numbers, API snapshot, downstream check script | Done (PR #6) |
| 3 | Member model in `structure/` (contributions, builders, assembler, span laws, span field nodes) | Done (PR #6) |
| 4 | Sections assembled through `structure/` (`section/assembly/plan.py`, `members.py`) | Done (PR #7) |
| 5 | Move to the target tree, one PR per top-level package; then `SectionSpec` and `SectionOptions` | **Done** — tagged `v0.5.0` |
| 6 | Seams with materials and meshing (thickness law protocol, named edges, forced stations, size fields, merged cap and skin option) | **Done.** Seams in this kernel; Gmsh and constitutive plies stay in consumers |
| 7 | Completeness (Phase 5 list below) | **Done** (7.1–7.10). Mid-surface patches and optional STEP (`[cad]`) shipped |
| 8 | Exactness (Phase 4 below). Tag 0.6.0 after mid-surface, plots extra, and alias removal | **Done.** Nose, insert polygons, linear loft scalars (`TableLaw`; PCHIP removed), mid-surface field from `shell_graph`, plots extra, and 0.5.x alias removal |
| 9 | Separate blade-model package | **Declined.** Consumers keep their own adapters |

Step 5 breakdown (Phase 2 complete; tagged `v0.5.0`):

| Slice | What | State |
|-------|------|-------|
| Phase 1a/1b | One shape model: CSG shapes are `SDFNode` subclasses, compile step removed | Done (PR #8) |
| Phase 1c | One polygon gateway: `section/boundary/polygons.py` (Shapely then; GEOS `PolygonOps` now) | Done (PR #9) |
| `sdf/`, `integration/`, `frames/` | Move to the target tree | Done (PR #10, #11, #12) |
| 5a | `props` → `properties`, `recipes` → `walls` | Done (PR #13) |
| 5b | `constructive` into `boundary` and `properties` (`export.py` → `boundary/walls.py` and `cores.py`; `labels` and `emitters` → `boundary/`; `properties.py` → `properties/polygons.py`) | Done (PR #14) |
| 5c | `placement` into `layout` | Done (PR #15) |
| 5d | `blade/` split `cut.py` (`shell_graph.py`) | Done (PR #16) |
| Phase 2 | `SectionSpec`, `SectionOptions`, `api` → `spec`; then tag 0.5.0 | **Done** (folder move + keyword aliases + wall-type registry). Tagged `v0.5.0` |

Decisions already made by the maintainer:

- Port the four-kind member model (wall partition, branch, region, transverse member),
  geometry only. Thickness is a law in metres; plies and layup libraries stay downstream.
- The parametric layout catalogue (`{X}{Y}-{Z}`, 40 layouts) stays the front door. Each
  id compiles to a list of contributions. Ribs are a second, spanwise axis of a concept,
  not a fourth letter in the code.
- Spar caps stay a separate stacked wall layer under the skin, so geometry is unchanged.
  A merged cap and skin laminate (one wall on the combined midline) becomes an opt-in
  section option in step 6.
- The public name for webs is `ShearWeb` (blade-anatomy term). `SpanWeb` was
  removed in 0.6.0.
- `sdf/` stays self-contained: nothing in it imports `section/` or `blade/`.
- `section/` stays one package, with subpackages in pipeline order (decided 6 Nov).
- Release tags are made by the maintainer only.
- This package is the source of truth for **blade geometry parameterisation**
  (airfoils, planform, layout catalogue, gauges in metres, members, span topology,
  `shell_graph`). It compiles published windIO (and parametric) inputs into
  `BladeGeometry` for suite benchmarks. Source YAML stays in the private
  reference-blade catalogue. Constitutive plies, loads, solvers and Gmsh stay in
  consumers. `src/` names no consumer.
- Plots live in `implicit_geometry_kernel.plotting` behind `[plots]`; root
  `plotting/` was removed in 0.6.0 (decision 1).
- The mid-surface query field is a view of `shell_graph`, not `(oml + cavity) / 2`
  as the spec (decision 2).
- Older demos under `examples/demos/` stay; port them onto `SectionSpec` /
  `SectionGeometry` / `BladeGeometry` / `shell_graph`. Do not delete (decision 3).
- Seams with materials and meshing live **in this kernel** (step 6). A Gmsh driver
  and ply ABD/PCOMP do not. No third blade-model package (decision 4; step 9 no).
- Public docs describe the three geometry products and roles. `src/` and `tests/`
  must not name consumer packages (`tests/test_consumer_names.py`). Scripts and
  `docs/FORK.md` may name checkouts.
- Kernel mesh is `build_shared_quadtree`. Optional extra `[mesh]` (`gmsh`) is
  leftover convenience, not advertised. `[cad]` / STEP patches shipped after
  step 6.
- 0.6.0 closes the deprecation window: mid-surface exactness, plots extra,
  alias removal. Maintainer tags it.

## How to work on a step

Every PR must show that it changes no numbers, except PRs in the exactness step, and
those state the size of each change in the PR body.

```bash
ruff check src tests scripts examples
python -m pytest -q                                # includes tests/test_snapshots.py
python scripts/snapshots.py                        # golden numbers and API, no update
python scripts/check_consumer.py PATH/TO/CONSUMER  # each downstream checkout
```

- **Golden numbers** (`tests/golden/*.json`, relative tolerance 1e-9): section
  properties, midline strips and cell loops for every catalogue layout at 0, 15 and 35
  degrees of twist, plus a lofted blade cut at five stations. Regenerate with
  `python scripts/snapshots.py --update` only in an exactness PR, and say why.
- **API snapshot** (`tests/golden/api.json`): new public names pass; a changed or removed
  signature fails. Regenerate only when the change is intended and in the changelog.
- **Downstream check**: run `scripts/check_consumer.py` against each package that pins
  the kernel. It installs this checkout and runs that package's tests.
- **One concern per PR.** Moves and renames go in their own PRs, apart from logic
  changes.
- **Deprecate, then remove.** An old import path stayed as an alias module for one minor
  release (pattern: `import <new> as _moved; sys.modules[__name__] = _moved`) and was
  listed in the changelog. Those 0.5.x aliases were removed in 0.6.0.
- Add an entry under `## Unreleased` in `CHANGELOG.md` for every PR.
- Never push to `main` directly.

## Target architecture

The kernel is a general-purpose library for blade and wing section geometry, and the
source of truth for that **parameterisation**. It names no consumer in `src/`: no
repository-specific code paths, options or defaults. A consumer pins a release tag,
imports the public API, and keeps adapters (mesher, ply stacks, solvers) on its own
side. The three products — integrated properties, wall midlines / `shell_graph`,
cell wall loops — are the shared geometry for beam FE, thin-wall section models,
constitutive thickness, and conformal shell FE.

Layers, with imports pointing down only:

```
contract  <-  sdf / integration / frames  <-  structure  <-  catalogue  <-  section  <-  blade  <-  extract
```

Extras (`plotting` `[plots]`, leftover `mesh` `[mesh]`, `blade.cad` `[cad]`,
`profiling.py`) are not layers. Geometry modules must not import them. The
in-kernel mesh is `integration/quadtree`. Consumers own shell mesher jobs;
`mesh/gmsh/` is unadvertised convenience.

- `catalogue/` holds WX-Y-Z-NR tokens; it imports only itself (and `contract`).
- `structure/` defines what members exist (wall partitions, branches, regions,
  transverse members) in geometric terms only.
- `section/` turns them into exact walls, midlines and properties at one station.
- `blade/` lofts stations and places transverse members along the span.
- Constitutive plies, Gmsh, solvers and loads stay out of the kernel. Step 6
  defines seams (`ThicknessLaw`, named edges, size fields) so consumers plug in.

### Robustness rules

| Area | Rule | Enforced by |
|------|------|-------------|
| Agnostic | No consumer names, paths or special cases in `src/` or tests. Docs describe products and roles | `tests/test_consumer_names.py` (src and tests) |
| Public API | Only names in the top-level `__all__` and `contract/` are public | API snapshot test |
| Versioning | Semantic versions; removals only in a minor release (pre-1.0) after one release of deprecation warnings; tags never move. Releases stay on 0.9.x; do not tag `v1.0.0` (maintainer, 11 Nov: the public API is not worked out) | Changelog; tags by the maintainer |
| Main branch | Every change through a PR with green CI (lint, tests on Python 3.11 to 3.13, examples, wheel build) | Branch protection on `main` |
| Numbers | Results reproducible between releases unless the changelog says otherwise | Golden-number test |
| Dependencies | Runtime: numpy, numba. Offset/overlay: ctypes ``libgeos_c`` (system/conda GEOS). Extras: scipy, matplotlib, cadquery, gmsh, pyyaml, contourpy | Core-install CI job |
| Typing | Ships `py.typed`; public API annotated | mypy on `src/` in CI |
| Data | No private reference-blade inputs in the repo; examples that need them skip without a local clone | Example 09 pattern |
| Submodule use | `pip install -e path/to/kernel`, no build step or generated files | Wheel-build job |

### Target tree

The end state. Every phase moves code toward it, and nothing is added outside it.

```
implicit_geometry_kernel/
  contract/        hand-off types (public)
  geos/            ctypes libgeos_c, WKB, kernel Polygon (C++ swap surface)
  sdf/             distance fields, self-contained
    nodes/           base, primitives, csg, bbox, tagged, span, lofted   (done)
    kernels/         numba kernels: primitives, polygon, grid, csg, transforms   (done)
    shapes/          airfoil and NACA fields   (done)
    solve/           GeometrySolver
  integration/     shoelace sums; quadtree/ (tree, soa, marching, moments, quadrature)   (done)
  frames/          section frames; tensors/ (inertia, moments, forces, thickness)   (done)
  structure/       contributions, builders, assembler, span laws   (done)
  catalogue/       WX-Y-Z-NR tokens (not inside section/layout)
  section/         in pipeline order:
    spec/            SectionSpec, SectionOptions, validation
    layout/          station layout, through-thickness anchoring   (done)
    walls/           one module per wall type, registered (skin, cap, web, core, insert)   (done)
    assembly/        plan, members, MultiCellSection
    boundary/        polygons, walls, cores, emitters, labels   (done)
    midlines/
    junctions/
    cells/
    properties/      including polygons.py from constructive   (done)
  blade/           loft façade, stations, cut, shell_graph, oml, concepts/, topology/
  extract/         marching (zero loops), curves (feature curves)
extras (not layers):
  plotting/        [plots]
  mesh/            leftover [mesh] convenience (not the kernel mesh)
  blade/cad.py     [cad] STEP writer
  profiling.py     optional timers, not on root __all__
```

Rules that keep it fixed:

- At most three levels: package/subpackage/module.py. Exception: a former oversized
  module may become a subpackage (`blade/topology/`, `sdf/kernels/`). No further nesting.
- Modules are 400 lines or fewer; packages hold 12 modules or fewer. A file that grows
  past the limit splits within its own folder.
- The import-boundary tests enforce the layer order.
- `tests/` and `docs/` mirror `src/` folder for folder.

Large modules and where they go:

| Module today | Split into |
|--------------|-----------|
| `integration/quadtree.py` | Done: `quadtree/tree.py`, `soa.py`, `marching.py`. Object-model `refine.py` / `refine_multi.py` dropped |
| `section/assembly/multicell.py` | webs, caps, skin and a short `section.py` (partly done by `plan.py` and `members.py`; ~261 lines) |
| `frames/tensors.py` | Done: `tensors/inertia.py`, `moments.py`, `forces.py`, `thickness.py` |
| `section/constructive/export.py` | Done: `boundary/walls.py`, `boundary/cores.py` |
| `integration/quadrature.py` | Done: `quadtree/moments.py`, `quadtree/quadrature.py` |
| `section/api/worker.py` | Done: `spec/build.py` (emit), `spec/prep.py` (cache key / sized assemble), `spec/worker.py` facade |
| `blade/cut.py` | Done: `cut.py` (SectionCut, walls, cells) and `shell_graph.py` (graph, partitions, topology diff) |
| `blade/loft.py` | Done: façade in `blade/loft/` (`loft.py`, `loft_tables`, `loft_cut`, `loft_fields`, `locked_loft`) |
| `blade/topology/map.py` | Done: `map.py` plus `map_walk` / `map_loft` |
| `blade/shell_graph.py` | Done: product plus `shell_graph_build` |
| `section/assembly/build_spar_caps.py` | Done: peel into `build_spar_caps_pads` |
| `plotting/plot_section.py` | Done: `plot_section_axes` / `plot_section` / `plot_section_hybrid` |

## Step 5 (continued): Phase 2 — one section input, no hidden settings

After 5c and 5d. This is still part of step 5, not a new status row.

- **2a.** `SectionSpec`, a frozen dataclass with one name per quantity: `layout`,
  `airfoil`, `twist_rad`, `gauges` (a `Gauges` object for skin SS/PS, cap and web,
  replacing both the thickness arguments and `zone_gauges_m`) and `webs` (x positions in
  metres, optional y anchors). It validates once in `__post_init__`; keyword alias
  rules live on `from_kwargs` / the helpers in `section_spec.py`.
- **2b.** `SectionOptions` for computational choices: outputs, props backend, quadtree
  settings (depth, max cells, band scale), sample counts, retain_mesh. The
  `SECTION_GEOMETRY_*` environment variables become defaults on this class (profiling
  may keep its variable; it does not change results).
- **2c.** `build_section_geometry(spec, options=SectionOptions())` is the only entry
  point. The keyword form becomes a thin wrapper that builds a spec and warns.
  `BladeGeometry.section_at` builds a `SectionSpec` internally.

## Step 5 (continued): Phase 2d — one place per wall type

**Done.** Station wall kinds dispatch on `SECTION_BUILDERS` (`on_shell_product`,
`polygon_by_surface`). `StiffenerWall` is a registered constructive kind. Label
helpers (`is_skin_label`) remain indexing, not a second kind registry.

Adding a member type means one module in `section/walls/` plus a registration,
with no edits elsewhere. Phase 2d **extends** the existing `structure` builders
and `SECTION_BUILDERS`. Do not add a second, parallel registry.

```python
class WallType(Protocol):
    kind: str
    def place(self, spec) -> list[WallInstance]                 # where it sits
    def field(self, inst) -> SDFNode                            # distance field
    def polygon(self, inst, boundary) -> Polygon                # exact solid
    def midline(self, inst, boundary) -> Polyline | None        # None for solids (cores)
    def junctions(self, inst, others) -> list[JunctionRule]     # how it meets other walls
```

- Each stage (assembly, boundary, midlines, junctions, cells, properties) loops over the
  registry and never tests names. Labels come from kind plus index, made in one place.
- Contract tests: a dummy stiffener type registered only in tests flows through to
  properties, midlines and cell loops without any edit in `src/`.

Member kinds (already in `structure/`):

| Kind | Lives in | What it is | Members |
|------|----------|------------|---------|
| Wall partition | section | extra thickness inside a footprint on the wall | skin, spar caps |
| Branch | section | own mid-surface meeting the wall; may have a free end and partial span | shear webs, stringers |
| Region | section | volume only, no mid-surface | cores, LE/TE inserts, bondlines |
| Transverse member | blade, at a station z | lies in the section plane | ribs, bulkheads, frames |

A rib is not a section wall type: it appears only in the cut at its own z. So it is a
second extension point at blade level (`TransverseMember`: name, z, region). A branch
with `extent=(z0, z1)` makes its ends topology events (existing `TopologyEvent`).

Concept model for the long term: `BladeConcept = SectionLayout + SpanwiseScheme`. The
scheme is a typed list of spanwise members (bulkhead, rib, frame, partial-span stringer)
with position and parameters. Short codes such as `2D-F` are named presets; the typed
spec is the source of truth.

## Phase 3: names a new reader can follow

Not a status-table step. 3d is already tracked as step 5.

- **3a.** Done: tracking codes in `src/` and `tests/` comments are plain sentences;
  history still lives in `docs/section/bug_log.md`.
- **3b.** Done: `Airfoil` is the shape; `AirfoilSDF` is its field via
  `Airfoil.to_sdf()`. NACA maths lives in `sdf.shapes.naca`.
- **3c.** Done: one name per concept; dual names removed in 0.6.0
  (`CellPolygon`, `interpolate`, `SectionGeometry`, `web_offsets`, `ShearWeb`).
- **3d.** The tree moves (step 5). Moves only, no logic.
- **3e.** Plots: `implicit_geometry_kernel.plotting` behind `[plots]`; root
  `plotting/` removed in 0.6.0.
- Also: Done: leftover downstream package names stripped from live `docs/`
  (fork history stays in `docs/FORK.md`).

## Step 6: seams with materials and meshing

The kernel defines small interfaces; the other side plugs in. The kernel imports
neither Gmsh nor a ply library. Adapters live in consuming tools, not in a
separate blade-model package.

| Seam | Kernel provides | Other side provides |
|------|-----------------|---------------------|
| Thickness | `ThicknessLaw` protocol: T at (z, s) in metres, plus footprint fields | a layup-based law built from a ply stack, with ply-drop tapers |
| Edges | partition and ply-drop edges as named curves at every station, plus forced stations (member ends, ribs) | the mesher conforms and maps names to physical groups |
| Size | distance, thickness and curvature fields | the mesher's size callback |

Also in this step: the merged cap and skin option (one wall laminate on the combined
midline instead of a stacked cap layer), opt-in, default off, golden numbers unchanged
when off.

## Step 7 (Phase 5): completeness

Each gap is closed with a worked example and a test. Order of need:

1. **Done.** Transverse members (ribs, bulkheads, frames) with rib faces in the
   shell graph (`TransverseMember`, `BladeGeometry.transverse`, `ShellGraph.faces`).
2. **Done.** Branches with free ends and partial span (stringers).
3. **Done.** CST airfoils, thickness-ratio families about the camber line, flatback/blunt TE with
   a set TE thickness, Selig and Lednicer readers, interpolation by t/c.
4. **Done.** Span laws linear (`TableLaw` / `np.interp`); PCHIP removed. Root cylinder to airfoil transition as a named law.
5. **Done.** Thickness fields over (z, s) with ply-drop ramps as geometry (law plus taper length).
6. **Done.** Named key points (LE, TE, cap edges, web feet) at every station.
7. **Done (OML).** Watertight OML surface mesh (STL/OBJ) from a quad span grid.
   Mid-surface patches per member and optional STEP behind `[cad]` shipped.
8. **Done.** Point queries: thickness, member and local frame at any 3-D point; distance
   to the OML; batched.
9. **Done.** Planform plus airfoil CSV reader and a windIO writer so a geometry round-trips.
10. **Done.** `validate()`: members leaving the span, too-thin walls, unclosed OML,
    with station and member named.

Out of scope on purpose: constitutive plies and stiffness; Gmsh; solvers; loads.
Geometry seams for thickness, named edges and size fields are step 6.

## Step 8 (Phase 4): exactness

Each change has a test that measures the error before and after against an analytic
reference. Golden numbers are regenerated in these PRs with the change stated.

| Gap | Today | Change | Check |
|-----|-------|--------|-------|
| Airfoil nose | surfaces looked up by x after a cubic resample, split at the smallest-x vertex (NACA 0015 upper surface 8e-5 c off at x = 0) | **Done.** split at the true LE (farthest from TE); NACA 4-digit thickness from the formula | error against analytic NACA thickness at x/c 0 to 0.01 below 1e-8 |
| Mid-surface field | `(outer + cavity) / 2`, only a true mid-surface where thickness is uniform | **Done.** signed SDF of the `shell_graph` skin loop | `|phi|` at loop-run interiors vs graph midpoints on all 40 layouts below 5e-4 |
| Quadtree remainder | LE/TE inserts still integrated on a depth-7 quadtree | **Done.** emit inserts as polygons; keep the quadtree as a test cross-check | insert area and moments match shoelace sums to 1e-12 |
| Loft scalars | chord, twist, pitch axis, offset linear between stations | **Done.** linear only (`TableLaw` / `np.interp`); PCHIP removed — sparse knots invented curvature. Airfoils stay linear blends | golden `blade.json` stations match linear tables |

Polygon distance (`sdf_polygon`, winding-number sign) is already exact and stays.

## Standing work (not new steps)

- Older `examples/demos/` stay on the gallery (decision 3). Teaching, section
  demos, contract 01, feature-field loft and box-beam assemble use `SectionSpec`
  / `TagNode`.
- Kernel station/CSG mesh is `build_shared_quadtree` (`integration/quadtree`).
  Do not advertise `[mesh]` (`gmsh`); it is leftover convenience. Consumers own
  shell mesher jobs and ply ABD/PCOMP.
- G10/G11 are **closed for this kernel** (thin-wall adapter and PCOMP stay in
  consumers). Remaining park: TE classifier on OML only, no `T` taper, no cores
  on `shell_graph` unless a shell consumer needs them
  ([MERGED_SHELL_GAPS](section/reports/MERGED_SHELL_GAPS.md)).
  `LoftedPolygonSDF` BVH (`accelerate` / `accel_at`) is shipped.
- Mapped loft freeze: loft webs follow intercept `DroppedWeb.index`; no
  per-blade branches in loft code ([blade bug log](blade/bug_log.md)).
- Efficiency is an **open scan**, not an infinite emit track
  ([PERFORMANCE_WORK_SUMMARY](section/reports/PERFORMANCE_WORK_SUMMARY.md)
  kernel scan 2026-11-10). Resume code only with a named timing; do not reopen
  Phase 15 Boolean shaving.

## Risks

- Silent number drift during moves: caught by the golden test.
- Downstream breaks: caught by `check_consumer.py` on each PR; consumers stay on tags.
- Large shape-model changes: done (Phase 1 landed in PR #8 without number changes).
