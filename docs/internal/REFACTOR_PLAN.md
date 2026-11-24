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

The repository follows the shared development standard (synced into
`docs/conventions/standard/`, values in `.github/repo.yml`) from 0.9.6; the
kernel was its pilot. Changelog headings now use `## [Unreleased]` /
`## [0.9.N] - YYYY-MM-DD`.

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
| K6 | Tip closure on every blade: closes the skin shell at z/L = 1 with the polygon mesher, skin partition; `nR` counts interior ribs only | Done (#54, 0.9.2) |
| K7 | Mesh product: skin/web quads, polygonal ribs and tip, arrays with partition and thickness, no Gmsh | Done (#55, 0.9.2) |
| K8 | Leading-edge and spanwise refinement controls | Done (#56, 0.9.2) |
| K9 | windIO read on the blade's own sections: arcs by OML arc length at every z (corner to corner), per-side cap windows, two-foot webs, web span windows (BG-036) whose ends are mesh stations, tip closure on single-cell tips; station airfoils that reproduce their points (blunt corner, round end, base listed as points) | Done (#58) |
| K5 | Trimmed top-level API (deprecate first), internal docs to `docs/internal/`; still 0.9.x | **Done** — `docs/internal/`, flat mesh shims removed, rare topology helpers demoted from top-level `__all__` with root `__getattr__` warnings (0.9.x). Stay 0.9.x, no 1.0 |
| 0.9.3 | Release prep: CI green on the split `mesh/` packages (lint, mypy, fan-budget test, example runner), doc links after the `docs/internal/` move, `__version__` 0.9.3. Consumers then pin `v0.9.3` instead of the `K9` tag | Done (#78, tagged v0.9.3) |
| 0.9.4 | Hardening: CI, docs site, module split, repository rules | Merged (#79); the `v0.9.4` tag sits on the 0.9.3 commit, so its code is 0.9.3 |
| Release fixes | `release.yml`: per-tag job concurrency, every commit in a push checked, hand run finishes a release whose tag exists (#89) | In review |
| 0.9.5 | Release of #79 and #81 (docs site polish, GitHub releases, CodeQL, typing, slow-test split); flat mesh shims removed in 0.9.6 | Merged (#81, 8f11101); tag v0.9.5 made by running the `release` workflow on 8f11101, then consumers move to it |
| 0.9.6 | SG-033 flapwise twist limit, SG-034 key-point LE/TE under twist, SG-018 re-check, SG-035 web feet built on the skin (no snaps), SG-037 OML offset, flat mesh shims removed | In progress |
| Manual | Docs as U, R, V, D manuals (issue #90): PR 1 restructure, R0 notation and bibliography, R1 frames; then R2–R3, R4–R6, R7–R10, V, U, D. Plan: issue #90 | PR 1 merged (#91); next R2–R3. CF-004 (moment transform twist sign) fixed in #93 |
| L1 | Decisions 0014 (smooth authored laws, PCHIP with a station-density check, t = f(s, z)), 0015 (partition mesher), 0016 (property-only sub-partitions), 0017 (analytic surface derivatives); issues #95, #96, #97, #98 | Done (#99) |
| L2 | `PchipLaw` and `span_interp="pchip"` (opt-in), leave-one-out station-density check that refuses or flags under-sampled inputs (#95) | Done (#100) |
| L3 | `ThicknessField(t, z)` per partition, per-node thickness on `ShellMesh` (#95) | Done (#101) |
| L4 | Tagged `Skeleton` view of the shell graph, normalised t per partition (#96) | Done (#101) |
| L5 | Generic `SubPartition` in `structure/`; caps re-expressed with goldens unchanged; then flanges, web foot, wedge (#97) | Part done (#103): `SubPartition` bands (t window, optional z window, thickness number, law or field, eccentricity vector), overlap check, `ShellGraph.cap_sub_partitions()`, `assign_sub_partitions` on the shell mesh with a straddle flag. Still to do: mesh lines at band edges (with L6), authored flange, web foot and wedge builders |
| L4b | Smooth section curves (C2 spline through given points, or CST) with the station-density check; analytic patch derivatives x_t, x_z, x_tt, x_tz, x_zz, normal, principal curvatures, signed-distance gradient and Hessian (#98) | Part done (#102): `PartitionPatch` x(t, z) with analytic x_t, x_z, x_tt, x_tz, x_zz, normal, principal curvatures, kind tags and Newton closest point. Still to do: smooth section curves (spline or CST) and the signed-distance Hessian |
| L6 | `shell_mesh(method="partition")` beside the ruled-strip mesher, compared in CI; terminations as transition bands; junction nodes by Newton on both stacked distance fields, with an on-both-surfaces invariant (#96, #98) | Part done (#104, #105): `partition_mesh` / `shell_mesh(method="partition")` for constant topology and a closed loop, skin and webs, nodes shared by skeleton identity, mesh lines at web feet and band edges, compared with the ruled-strip mesh in tests; rib and tip-closure triangle fills on the station nodes (#105); `method="auto"` falls back to the ruled-strip mesher where the partition mesher refuses (#106). Skin node spacing follows section curvature (sagitta rule, graded) and span stations follow the skin edge (#114). Still to do (later, maintainer's choice 2026-11-24): transition bands at terminations, branches, Newton junction nodes between stations |
| L6b | Nodal and element curvature (k1, k2, theta) and kind tags on `ShellMesh`; kind by construction with exact zeros, computed curvature checked against the tag (#98) | Done (#102, #105): element curvature and kind (`element_curvature`); nodal values per partition (`node_curvature`), evaluated at each node's own (t, z) |
| L7 | Flagged release: PCHIP default for authored blades, partition mesher default, goldens re-baselined (#95, #96) | Done (#113, maintainer's approval 2026-11-24): `span_interp="auto"` (PCHIP per law where the leave-one-out check passes, linear otherwise; `from_section` linear) and `shell_mesh(method="auto")` (partition mesher when sized by `h` alone, ruled-strip fallback) are the defaults; golden numbers unchanged, so no re-baseline was needed. Ruled-strip mesher kept for one more release |
| L8 | Alongside items from the 2026-11-22 review: mesh provenance (`ShellMesh.digest`, `provenance`), shell-mesh benchmark and binary VTU, junction mode B (one foot plus angle), strict invariant gate and property fuzzer | Done (#107, #109, #110, #111): provenance; binary VTU and `scripts/bench_shell_mesh.py`; junction mode B (`web_tilt_anchor`); strict gate (`check_shell_mesh`, `shell_mesh(strict=True)`) and Hypothesis fuzzer |

#### Consumer adapter gaps (to do)

Maintainer rule (12 Nov): a consumer hands all of its geometry to the kernel. Anything
a consumer's adapter still has to compute, convert, approximate or refuse is a kernel
to-do, listed here until the adapter can drop it. The list comes from building a kernel
blade out of a parametric shell-model consumer's inputs (law-defined planform, airfoil
family with a root-circle blend, webs, caps, ribs, stringers) and out of windIO files.

| Id | What the adapter does today | What the kernel should do | State |
|----|-----------------------------|---------------------------|-------|
| G1 | Maps windIO `nd_arc` itself, because the reader used `x/c = 1 - 2 min(s, 1 - s)` (webs 25–290 mm off on a 15 MW blade) | Place arcs on the section the blade has at that span (`blade/windio_arc.py`) | Done (K9, #58) |
| G2 | Refuses suction and pressure caps with different chord bands | Per-side cap windows (`CapPlacementRule.from_side_windows`, `cap_windows` returning a mapping) | Done (K9, #58) |
| G3 | Cannot pass a web through its two attachment points (one chord position per web) | Two-foot webs: `web_positions` pairs of (suction, pressure) fractions | Done (K9, #58) |
| G4 | Works round partial-span webs kept root to tip (BG-036) | Web span window from the `nd_arc` grid; layout steps down where no web is live | Done (K9, #58) |
| G5 | — (found on the way) | A single-cell tip at a coarse `h` lost its tip closure (4-node ring); the loop now closes on at least five nodes | Done (K9, #58) |
| G5b | Could not get within 0.1 mm at a trailing edge however many stations it added | Station airfoils pass through their given points; round trailing ends and bases listed as points are kept (BG-042, BG-043) | Done (K9, #58) |
| G5c | Saw webs start and end at the nearest mesh station | Partial-span web ends are forced stations (BG-045) | Done (K9, #58) |
| G6 | Bisects stations until the loft is within 0.1 mm of its own analytic surface (81 stations become a few hundred) | Take the outer shape as laws or a callable `section(z)` and refine stations to a stated tolerance | Done (`BladeGeometry.from_section`; linear samples, no second morph) |
| G7 | Adds quarter-point probes because the airfoil morph changes rule (smoothstep or linear) when one station is at least 60 % thick and the next is not | One documented morph rule, or the rule as an explicit option | Done (`airfoil_morph`: auto, smoothstep, thickness, linear) |
| G8 | Converts every chord fraction: its fractions are measured from the pitch axis on the chord law, the kernel's on each station's own x extent (they differ across a root blend) | Accept positions in metres from the pitch axis, or a reference chord law | Done (`web_x_from_pitch`, `cap_x_from_pitch`, `Stringer.position_m`, `chord_fraction_from_pitch`) |
| G9 | Flips the twist sign | State the twist convention in `Station` and offer the other sense | Done (`twist_sense`: `nose_up` default, `nose_down`) |
| G10 | Leaves `merge_cap_and_skin` off: with it on the caps vanish from the shell graph and one cap run meshes a degenerate triangle | Fix; add to the bug log | Done (BG-040: keep the cap partitions on the merged cut) |
| G11 | Records a cap-edge difference (up to 17 mm on a root circle): the kernel cuts the edge at `x` on the cap mid-plane, the consumer at `x` on its node surface | Choice of the surface the edge is measured on | Done (`cap_edge_reference`: `midplane` default, `outer`) |
| G12 | Refuses tilted webs | A tilt angle on a web (two-foot webs already carry the geometry) | Done (`web_tilt_deg` on a chord-fraction web) |
| G13 | Refuses webs with different layups; one cap thickness for both sides | Thickness per web and per cap side | Done (`spar_cap_thickness=(suction, pressure)`, `shear_web_thickness` one value per web) |
| G14 | Refuses wall nodes on the outer surface (`wall_reference = oml`) | A reference-surface option for the shell graph (mid-surface, OML) | Done (`wall_reference`: `mid` default, `oml`) |
| G15 | Refuses stringers whose chord position or height varies along the span | Stringer position and height as span laws | Done (`Stringer.position_law`, `Stringer.depth_law`) |
| G16 | Finds windIO spar caps by layer name (`spar`, `cap`); a file that names layers by keypoint (10 MW) gets none | Identify caps by structure (the band between the web feet, uniaxial material) or take an explicit list | Done (`cap_id`: name rule default, `between_webs`, or a name list; BG-041) |
| G17 | Sums ply thicknesses itself for skin, cap and web | Thickness per partition from the stack (decision D4) | Done (`StackThickness.from_layup`, `StackThickness.from_windio`) |
| G18 | — | Webs that cross between stations (BG-037) | Done (web anchors follow the placed arcs) |
| G19 | Splits rib and tip-closure polygon cells into triangles for a solver element that takes triangles. A fan split at a reflex corner put a triangle outside its cell, over the next one (0.09 % extra rib and tip area at a production size; nothing raised because each triangle was rewound to positive area) | Ship a valid split of each polygon cell (through reflex corners, every piece inside its cell, a piece outside raises), so the consumer takes cells or triangles from the kernel | Done (0.9.3 meshes rib and tip cells as triangles, unstructured, with boundary recovery; `tests/test_rib_triangles_cover_cells.py` checks no fold and no overlap per station plane) |

Not kernel work: the trailing-edge bond and core fill are materials inside the section,
and plies, loads and solvers stay in consumers.

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
- Release tags are made by the `release` workflow when a release PR merges
  (decision 0013, 22 Nov; it replaced "maintainer only").
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
| Main branch | Every change through a PR with green CI (lint, doc links, tests on Python 3.11 to 3.13, examples, wheel build, consumer suites) | Branch protection on `main`, no bypass |
| Releases | A `v*` tag equals `v` + `__version__`; wheels go to PyPI only when publishing is switched on | `release.yml` |
| Docs site | Manual (the Markdown tree) plus an API reference from each package's `__all__`; built on every PR with warnings as errors, published to GitHub Pages from `main` | `docs.yml` |
| Contributor files | CITATION.cff, CONTRIBUTING, SECURITY, issue and PR templates, Dependabot, pre-commit, API stability policy; rules in AGENTS.md, scope in SCOPE.md, decisions in docs/decisions/, CODEOWNERS | CI coverage report with a floor |
| Module size | `mesh/shell/shell_mesh.py` (2625 lines) and `mesh/tables/tables.py` (1238) split by job in 0.9.4, no behaviour change; old module names still resolve every helper | mypy `check_untyped_defs` on `mesh/` |
| Figures | PNGs rewritten only past 0.5 % changed pixels; mesh exports untracked; cross-machine font differences still exceed that, so regenerated PNGs are committed only when the example's subject changed | `figure_io.write_png` |
| Numbers | Results reproducible between releases unless the changelog says otherwise | Golden-number test |
| Dependencies | Runtime: numpy, numba. Offset/overlay: ctypes ``libgeos_c`` (system/conda GEOS). Extras: scipy, matplotlib, cadquery, gmsh, pyyaml, contourpy | Core-install CI job |
| Typing | Ships `py.typed`; public API annotated; strict (`disallow_untyped_defs`) on structure, contract, frames, catalogue and extract, roadmap in `docs/conventions/TYPING.md` | mypy on `src/` in CI |
| Docstrings | Public functions and classes documented; floor 69 % with plotting excluded, raised as coverage grows | `interrogate` in the lint job |
| Tests | Tests of 4 s or more are marked `slow` in `tests/conftest.py`; PRs run the full suite on Python 3.13 and fast tests elsewhere | `tests.yml` |
| CI installs | Dev and docs extras resolve from `requirements/ci.txt`; the core-install job stays unpinned | `tests.yml`, `docs.yml` |
| Releases page | A tag makes a GitHub release with its changelog section and the wheel | `release.yml` `github-release` |
| Security scan | CodeQL on Python and workflows, weekly and on PRs | `codeql.yml` |
| Traceability | Issue first, PR says `Closes #N`, changelog lines end with the PR number | `docs/conventions/standard/workflow.md`, PR template |
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
| Loft scalars | chord, twist, pitch axis, offset linear between stations | **Done.** linear only (`TableLaw` / `np.interp`); PCHIP removed — sparse knots invented curvature. Airfoils stay linear blends. Superseded for authored laws by [0014](../decisions/0014-smooth-authored-span-laws.md) (L2, L7) | golden `blade.json` stations match linear tables |

Polygon distance (`sdf_polygon`, winding-number sign) is already exact and stays.

## Standing work (not new steps)

- Older `examples/demos/` stay on the gallery (decision 3). Teaching, section
  demos, contract 01, feature-field loft and box-beam assemble use `SectionSpec`
  / `TagNode`.
- Kernel station/CSG mesh is `build_shared_quadtree` (`integration/quadtree`).
  Do not advertise `[mesh]` (`gmsh`); it is leftover convenience. Consumers own
  shell mesher jobs and ply ABD/PCOMP.
- G10/G11 are **closed for this kernel** (thin-wall adapter and PCOMP stay in
  consumers). The TE classifier park is closed (SG-032). `StackThickness`
  ramps at cap edges; cores sit on default `shell_graph.solids`
  ([MERGED_SHELL_GAPS](section/reports/MERGED_SHELL_GAPS.md)).
  `LoftedPolygonSDF` BVH (`accelerate` / `accel_at`) is shipped.
- Mapped loft freeze: loft webs follow intercept `DroppedWeb.index`; no
  per-blade branches in loft code ([blade bug log](../blade/bug_log.md)).
- Efficiency is an **open scan**, not an infinite emit track
  ([PERFORMANCE_WORK_SUMMARY](section/reports/PERFORMANCE_WORK_SUMMARY.md)
  kernel scan 2026-11-10). Resume code only with a named timing; do not reopen
  Phase 15 Boolean shaving.

## Risks

- Silent number drift during moves: caught by the golden test.
- Downstream breaks: caught by the `consumers` CI job (`check_consumer.py` on each repository listed in the `CONSUMER_REPOS` variable) on each PR; consumers stay on tags.
- Large shape-model changes: done (Phase 1 landed in PR #8 without number changes).
