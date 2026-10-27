# api — public entrypoints

↑ [`section_geometry`](../README.md)

Public surface: construct `SectionGeometryRunner`, optionally `resolve()` the
layout, `execute()` the build, read `SectionGeometryResult`. Validates inputs,
sequences work, and packages the product; domain maths lives in the other
modules.

Lifecycle: construct stores arguments; `resolve()` → `StationLayoutResolution`
(no solids); `execute()` builds gated outputs; `get_results()` returns the
frozen product (or raises, e.g. `InfeasibleSectionGeometry`). Sizing: resolve
once at min skin, pass `resolution=` on later trials.

Spine context: [`CASCADE_GLOSSARY`](../../../docs/CASCADE_GLOSSARY.md) ·
[`COEXISTENCE_IMPLICIT_SHELL`](../../../docs/COEXISTENCE_IMPLICIT_SHELL.md).

## Table of contents

1. [Which entrypoint](#which-entrypoint)
2. [Organisation](#organisation)
3. [API](#api)
4. [Caller parameters](#caller-parameters)
5. [Errors and feasibility](#errors-and-feasibility)
6. [Usage in the package](#usage-in-the-package)
7. [See also](#see-also)

---

## Which entrypoint

| Call | When |
|------|------|
| `SectionGeometryRunner` | Public call path |
| `build_section_geometry` | Free-function / tests |
| `build_section_with_thicknesses` | Assemble solids only (any key → `MultiCellSection`) |
| `build_section_view` | Demos / properties report (same MultiCell path; long thickness kwargs) |

## Organisation

| File | Role |
|------|------|
| `runner.py` | `SectionGeometryRunner` |
| `worker.py` | `build_section_geometry`, `resolve_station_layout`, result dataclasses |
| `factory.py` | `build_section_view` and related helpers |

## API

| Name | File | Role |
|------|------|------|
| `SectionGeometryRunner` | `runner.py` | Preferred lifecycle |
| `build_section_geometry` | `worker.py` | Free-function build |
| `resolve_station_layout` | `worker.py` | Drop webs that do not fit |
| `StationLayoutResolution` | `worker.py` | Effective layout + survivors + drops |
| `SectionGeometryResult` | `worker.py` | Return value |
| `section_geometry_prep_key` | `worker.py` | Cache key for a build |
| `build_section_view` | `factory.py` | Convenience factory for demos |

Product field tables: [`../README.md` § Results](../README.md#results).

## Caller parameters

Important keyword arguments on `SectionGeometryRunner` /
`build_section_geometry`:

| Parameter | Role |
|-----------|------|
| `layout` | Nominal `SystemLayoutSpec` (from `resolve_system_type`) |
| `airfoil` | Station `AirfoilSDF` |
| `skin_thickness`, `spar_cap_thickness`, `shear_web_thickness` | Geometric wall thicknesses [m] (material-bind agnostic). Short `t_skin` / `t_spar` / `t_web` are **not** accepted. Optional resolve aliases: `cap_height` ↔ `spar_cap_thickness`, `web_thickness` ↔ `shear_web_thickness` |
| `twist_rad` | Geometric twist [rad] for emit |
| `outputs` | `"properties"` / `"midlines"` / `"cell_wall_loops"` |
| `web_x_positions_m` | Absolute web x offsets from LE [m], fixed across span. If omitted, positions are `frac × local_chord` and webs are never dropped for taper |
| `resolution` | Precomputed `StationLayoutResolution` — skip re-resolve |
| `min_cavity_depth_m` | Drop a web whose local cavity is too shallow for min-gauge caps |
| `airfoil_for_spec` | Aerofoil used to size the quadtree domain; pass the **twisted** outline so the AABB does not clip rotated tips |
| `props_backend` | `"hybrid"` / `"quadtree"` / `"grid"` (default from env) |
| `quadtree_spec` | Optional QT refinement settings |
| `retain_mesh` | **G-PROPS debug only** — keep hybrid/quadtree occupancy mesh on the result. Not a `section_mesh` product |

Layouts `"0A"` / `"0B"`: valid zero-web hosts; shell-mesh handoff accepts
skin ± optional caps (no web junctions) — see
[`ADR_ZERO_WEB_SHELL_MESH.md`](../docs/reference/ADR_ZERO_WEB_SHELL_MESH.md).

## Errors and feasibility

| Situation | Behaviour |
|-----------|-----------|
| Some webs do not fit the local cavity | Those webs are **dropped**; `StationLayoutResolution.dropped` lists them; the effective layout may step down the catalogue |
| **No** web can be anchored (and the layout needs webs) | Raises `InfeasibleSectionGeometry` (≠ “too thick”) |
| Unknown `outputs` flag | Raises `ValueError` |

Thicknesses are clamped to package limits before assembly; clamping alone does
not raise.

Environment knobs: `SECTION_GEOMETRY_PROPS_BACKEND=hybrid|quadtree|grid`;
`SECTION_GEOMETRY_PROFILE=1` fills `result.stage_timings`.

## Usage in the package

```mermaid
flowchart LR
  Runner[SectionGeometryRunner] --> Resolve[resolve]
  Resolve --> Eff[StationLayoutResolution]
  Runner --> Exec[execute]
  Exec --> Build[build_section_geometry]
  Build --> Result[SectionGeometryResult]
```

## See also

- [`../layout/`](../layout/README.md) — layout resolve rules
- [`../types/`](../types/README.md) — public aliases
- [`../docs/guides/caller_patterns.md`](../docs/guides/caller_patterns.md) — sizing / spanwise patterns
