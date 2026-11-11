# Migrating from 0.2.x to 0.9

Code written against 0.2.x stops at the first of these changes. Each row gives
the 0.2.x call and its 0.9 form. Geometry is unchanged by the migration: on a
NACA 0024 two-web blade, web lengths, cap length and cell areas agree with
0.2.0 to 2e-4 relative (`tests/test_compat_0_2.py`).

## Install

| 0.2.x | 0.9 |
|---|---|
| Shapely 2 pulled in by pip | Install GEOS first: `apt install libgeos-dev`, `brew install geos` or `conda install geos`. The kernel loads `libgeos_c` through ctypes; set `GEOS_LIBRARY_PATH` if it is not on the loader path |

## Layout codes

A layout code names the section family, the web anchoring and the rib count.
The short form is rejected with `KeyError`.

| 0.2.x | 0.9 |
|---|---|
| `"2D-F"` | `"2D-F-LC-0R"` |
| `"0A"`, `"3B-CN"`, ... | `"0A-LC-0R"`, `"3B-CN-LC-0R"`, ... |

Anchoring is `LC` (local chord fraction), `LE` (distance from the leading
edge) or `PA` (blade-x from the pitch axis). Web positions given as chord
fractions match 0.2.x with `LC`. `0R` means no ribs; `10R` adds ten ribs.

## BladeGeometry

| 0.2.x | 0.9 |
|---|---|
| `BladeGeometry(stations, "2D-F", t_skin, t_cap, t_web, ...)` | `BladeGeometry(stations, "2D-F-LC-0R", t_skin, t_cap, t_web, ..., pivot_fraction=0.333)` |
| no `pivot_fraction` (fell back to 0.333) | `pivot_fraction` is required: a chord fraction or a span table. Use the stations' pitch axis where you set one |

## Imports

| 0.2.x | 0.9 |
|---|---|
| `implicit_geometry_kernel.section_geometry.layout.system_layout.SYSTEM_TYPE_KEYS` | `implicit_geometry_kernel.SYSTEM_TYPE_KEYS` |
| `implicit_geometry_kernel.blade_geometry` | `implicit_geometry_kernel.blade` |

## Section level

| 0.2.x | 0.9 |
|---|---|
| `SectionGeometry(layout=..., airfoil=..., ...)` | `SectionGeometry(SectionSpec.from_kwargs(...), SectionOptions.from_kwargs(...))` |
| `apply_shapely` | `apply_clip` |
| `cut.mid_surface()` put a cap + skin run on the stack's mid-surface | The run stays on the skin mid-surface. Give the stacked laminate an offset of half the cap thickness, inward (`cut.shell_graph().offset_regions()`) |

Changes made after 0.2 and removed again before 0.9 (`span_kind`, `SpanWeb`,
`ExprSDF`, `WallRegionGraph.cap_mode`) are listed in the changelog.

## Checking a consumer

Run the consumer's tests against a kernel checkout:

```bash
python scripts/check_consumer.py PATH/TO/CONSUMER
```
