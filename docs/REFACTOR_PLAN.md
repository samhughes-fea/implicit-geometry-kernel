# Refactor and exactness plan

This file is the working plan for restructuring the kernel and closing its exactness
gaps. It is written so that anyone (a person or an automated workflow) can pick the work
up from the **Status** section without other context. Update the status table in the
same PR that finishes a step.

## Status (resume here)

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
| — | Phase 3a–c, 3e | Naming and plots; 3e is blocked on open decision 1 |

| Step | What | State |
|------|------|-------|
| 1 | Stabilise: lint fix, naming fix, data clean-up | Done (PR #4, #5). Branch protection, consumer-name grep and mypy on `src/` added. `v0.4.1` was skipped; tagged `v0.5.0` |
| 2 | Safety net: golden numbers, API snapshot, downstream check script | Done (PR #6) |
| 3 | Member model in `structure/` (contributions, builders, assembler, span laws, span field nodes) | Done (PR #6) |
| 4 | Sections assembled through `structure/` (`section/assembly/plan.py`, `members.py`) | Done (PR #7) |
| 5 | Move to the target tree, one PR per top-level package; then `SectionSpec` and `SectionOptions` | **Done** — tagged `v0.5.0` |
| 6 | Seams with materials and meshing (thickness law protocol, named edges, forced stations, size fields, merged cap and skin option) | Not started. Coupling location is open decision 4 |
| 7 | Completeness (Phase 5 list below) | Not started |
| 8 | Exactness (Phase 4 below). Tag 0.6.0 after | Not started. Mid-surface row is blocked on open decision 2 |
| 9 | Separate blade-model package (only if the maintainer chooses it) | Not decided (open decision 4) |

Step 5 breakdown (Phase 2 complete; tagged `v0.5.0`):

| Slice | What | State |
|-------|------|-------|
| Phase 1a/1b | One shape model: CSG shapes are `SDFNode` subclasses, compile step removed | Done (PR #8) |
| Phase 1c | One Shapely gateway: `section/boundary/polygons.py` | Done (PR #9) |
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
- The public name for webs is `ShearWeb` (blade-anatomy term). `SpanWeb` is an alias
  until 0.6.0.
- `sdf/` stays self-contained: nothing in it imports `section/` or `blade/`.
- `section/` stays one package, with subpackages in pipeline order (decided 6 Nov).
- Release tags are made by the maintainer only.

Open decisions (ask the maintainer; do not assume; do not start the blocked work):

1. Plots package: rename it to `implicit_geometry_kernel_plots`, or move it back inside
   the kernel behind the `[plots]` extra? Blocks Phase 3e.
2. Phase 4 mid-surface: rebuild it from the polygons, or keep the field and add a drift
   check? Blocks the mid-surface row of step 8.
3. Should any of the 41 older demos in `examples/demos/section` be retired rather than
   moved? Blocks relocating those demos.
4. Where the coupling with materials and Gmsh lives (recommended: a separate blade-model
   package that imports the kernel; the kernel stays free of Gmsh). Blocks step 9 and
   the packaging of step 6.

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
- **Deprecate, then remove.** An old import path stays as an alias module for one minor
  release (pattern: `import <new> as _moved; sys.modules[__name__] = _moved`) and is
  listed in the changelog. Removal happens in 0.6.0.
- Add an entry under `## Unreleased` in `CHANGELOG.md` for every PR.
- Never push to `main` directly.

## Target architecture

The kernel is a general-purpose library for blade and wing section geometry. It names no
consumer: no repository-specific code paths, options, defaults or docs. A consumer pins a
release tag (as a git submodule or with pip), imports only the public API, and keeps its
own adapters on its own side.

Layers, with imports pointing down only:

```
contract  <-  sdf, integration, frames  <-  structure  <-  section  <-  blade  <-  extract
```

- `structure/` defines what members exist (wall partitions, branches, regions,
  transverse members) in geometric terms only.
- `section/` turns them into exact walls, midlines and properties at one station.
- `blade/` lofts stations and places transverse members along the span.
- Materials, plies, meshing, solvers and loads stay out of the kernel.

### Robustness rules

| Area | Rule | Enforced by |
|------|------|-------------|
| Agnostic | No consumer names, paths or special cases in `src/`, tests or docs. Generic extension points only | `tests/test_consumer_names.py` (src and tests) |
| Public API | Only names in the top-level `__all__` and `contract/` are public | API snapshot test |
| Versioning | Semantic versions; removals only in a minor release (pre-1.0) after one release of deprecation warnings; tags never move | Changelog; tags by the maintainer |
| Main branch | Every change through a PR with green CI (lint, tests on Python 3.11 to 3.13, examples, wheel build) | Branch protection on `main` |
| Numbers | Results reproducible between releases unless the changelog says otherwise | Golden-number test |
| Dependencies | Runtime: numpy, shapely, numba. Extras: scipy, matplotlib, pyyaml, gmsh | Core-install CI job |
| Typing | Ships `py.typed`; public API annotated | mypy on `src/` in CI |
| Data | No private reference-blade inputs in the repo; examples that need them skip without a local clone | Example 09 pattern |
| Submodule use | `pip install -e path/to/kernel`, no build step or generated files | Wheel-build job |

### Target tree

The end state. Every phase moves code toward it, and nothing is added outside it.

```
implicit_geometry_kernel/
  contract/        hand-off types (public)
  sdf/             distance fields, self-contained
    nodes/           base, primitives, csg, bbox, tagged, span, lofted   (done)
    kernels/         numba kernels: primitives, polygon, grid, csg, transforms   (done)
    shapes/          airfoil and NACA fields   (done)
  integration/     shoelace sums; quadtree/ (tree, refine, refine_multi, soa, marching,
                   moments, quadrature)   (done)
  frames/          section frames; tensors/ (inertia, moments, forces, thickness)   (done)
  structure/       contributions, builders, assembler, span laws   (done)
  section/         in pipeline order:
    spec/            SectionSpec, SectionOptions, validation (from api/worker, kwarg_aliases)
    layout/          catalogue, layout resolution, through-thickness anchoring   (done)
    walls/           one module per wall type, registered (skin, cap, web, core, insert)   (done)
    assembly/        plan, members, MultiCellSection
    boundary/        polygons, walls, cores, emitters, labels   (done)
    midlines/
    junctions/
    cells/
    properties/      including polygons.py from constructive   (done)
  blade/           loft, stations, cut, shell_graph, transverse members
  extract/         curves and surfaces
```

Rules that keep it fixed:

- At most three levels: package/subpackage/module.py.
- Modules are 400 lines or fewer; packages hold 12 modules or fewer. A file that grows
  past the limit splits within its own folder.
- The import-boundary tests enforce the layer order.
- `tests/` and `docs/` mirror `src/` folder for folder.

Large modules and where they go:

| Module today | Split into |
|--------------|-----------|
| `integration/quadtree.py` | Done: `quadtree/tree.py`, `refine.py`, `refine_multi.py`, `soa.py`, `marching.py`. Still to do (a logic PR, not a move; keep it out of step 5): keep the array-based path (`soa.py`) only and drop the object-model duplicates (`refine.py`, `refine_multi.py`) once nothing uses them |
| `section/assembly/multicell.py` | webs, caps, skin and a short `section.py` (partly done by `plan.py` and `members.py`; ~261 lines) |
| `frames/tensors.py` | Done: `tensors/inertia.py`, `moments.py`, `forces.py`, `thickness.py` |
| `section/constructive/export.py` | Done: `boundary/walls.py`, `boundary/cores.py` |
| `integration/quadrature.py` | Done: `quadtree/moments.py`, `quadtree/quadrature.py` |
| `section/api/worker.py` | Done: `spec/build.py` (emit), `spec/prep.py` (cache key / sized assemble), `spec/worker.py` facade |
| `blade/cut.py` | Done: `cut.py` (SectionCut, walls, cells) and `shell_graph.py` (graph, partitions, topology diff) |

## Step 5 (continued): Phase 2 — one section input, no hidden settings

After 5c and 5d. This is still part of step 5, not a new status row.

- **2a.** `SectionSpec`, a frozen dataclass with one name per quantity: `layout`,
  `airfoil`, `twist_rad`, `gauges` (a `Gauges` object for skin SS/PS, cap and web,
  replacing both the thickness arguments and `zone_gauges_m`) and `webs` (x positions in
  metres, optional y anchors). It validates once in `__post_init__`; keyword alias
  rules live on `from_kwargs` / the helpers in `section_spec.py`. `kwarg_aliases.py`
  re-exports them until 0.6.0.
- **2b.** `SectionOptions` for computational choices: outputs, props backend, quadtree
  settings (depth, max cells, band scale), sample counts, retain_mesh. The
  `SECTION_GEOMETRY_*` environment variables become defaults on this class (profiling
  may keep its variable; it does not change results).
- **2c.** `build_section_geometry(spec, options=SectionOptions())` is the only entry
  point. The keyword form becomes a thin wrapper that builds a spec and warns.
  `BladeGeometry.section_at` builds a `SectionSpec` internally.

## Step 5 (continued): Phase 2d — one place per wall type

Adding a member type should mean one module in `section/walls/` plus a registration,
with no edits elsewhere. Today kinds are still found by string prefixes (`is_skin_label`,
`startswith("shear_web_")`, `startswith("spar_cap")`) across `boundary/labels.py`,
`boundary/walls.py`, `boundary/cores.py`, `boundary/emitters.py`, `assembly/multicell.py`,
`midlines/export.py`, `properties/export_section_json.py` and `blade/cut.py` (old
`constructive/` and `props/` paths remain aliases until 0.6.0).

Phase 2d **extends** the existing `structure` builders and `SECTION_BUILDERS`. Do not
add a second, parallel registry.

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

- **3a.** Replace internal tracking codes in code comments (G-POLY, G-MID, SG-009,
  BL-028, "Wave 0") with plain sentences; link `docs/section/bug_log.md` where history
  matters.
- **3b.** `Airfoil` is the shape; `AirfoilSDF` is its field, made through
  `Airfoil.to_sdf()`. NACA maths lives in one place.
- **3c.** One name per concept; deprecate the other: CellPolygon/CellLoop,
  interpolate/station_at, SectionGeometry/SectionGeometryRunner,
  web_offsets/web_x_positions_m, SpanWeb/ShearWeb. Remove in 0.6.0.
- **3d.** The tree moves (step 5). Moves only, no logic.
- **3e.** Plots package (open decision 1). Blocked until the maintainer chooses.
- Also: remove leftover mentions of downstream package names from `docs/`.

## Step 6: seams with materials and meshing

The kernel defines small interfaces; the other side plugs in. The kernel imports
neither side. Where the other side lives is open decision 4.

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

1. Transverse members (ribs, bulkheads, frames) with rib faces in the shell graph.
2. Branches with free ends and partial span (stringers).
3. CST airfoils, thickness-ratio families about the camber line, flatback/blunt TE with
   a set TE thickness, Selig and Lednicer readers, interpolation by t/c.
4. Span laws with PCHIP for planform; root cylinder to airfoil transition as a named law.
5. Thickness fields over (z, s) with ply-drop ramps as geometry (law plus taper length).
6. Named key points (LE, TE, cap edges, web feet) at every station.
7. Watertight OML surface mesh (STL/OBJ), quad surface grid, mid-surface patches per
   member with shared edges; optional STEP behind a `[cad]` extra.
8. Point queries: thickness, member and local frame at any 3-D point; distance to any
   named member; batched.
9. Planform plus airfoil CSV reader and a windIO writer so a geometry round-trips.
10. `validate()`: overlaps, members leaving the OML, too-thin walls, self-intersections,
    with station and member named.

Out of scope on purpose: materials, plies and stiffness; meshing; solvers; loads.

## Step 8 (Phase 4): exactness

Each change has a test that measures the error before and after against an analytic
reference. Golden numbers are regenerated in these PRs with the change stated.

| Gap | Today | Change | Check |
|-----|-------|--------|-------|
| Airfoil nose | surfaces looked up by x after a cubic resample, split at the smallest-x vertex (NACA 0015 upper surface 8e-5 c off at x = 0) | split at the spline's true LE (farthest point from TE), look up by arc length, evaluate NACA from the formula | error against analytic NACA thickness at x/c 0 to 0.01 below 1e-8 |
| Mid-surface field | `(outer + cavity) / 2`, only a true mid-surface where thickness is uniform | build the midline from the exact wall polygons and use the field only for queries, or add a gradient drift check (open decision 2) | midline against polygon midline on all 40 layouts, worst gap reported |
| Quadtree remainder | LE/TE inserts still integrated on a depth-7 quadtree | emit inserts as polygons; keep the quadtree as a test cross-check | insert area and moments match shoelace sums to 1e-12 |
| Loft scalars | chord, twist, pitch axis, offset linear between stations | optional PCHIP (monotone, no overshoot); airfoils stay linear blends | reference planforms match their tables with smooth first derivatives |

Polygon distance (`sdf_polygon`, winding-number sign) is already exact and stays.

## Risks

- Silent number drift during moves: caught by the golden test.
- Downstream breaks: caught by `check_consumer.py` on each PR; consumers stay on tags.
- Large shape-model changes: done (Phase 1 landed in PR #8 without number changes).
