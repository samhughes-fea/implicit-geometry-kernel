# Graphify EXTRACTED — no next kernel step

**Audience: internal.** Maintainer notes; not part of the library manual.

Dated 2026-11-10. Answers only. No `src/` change in this note.

Scan: `--update` on `src/implicit_geometry_kernel`. **3064 nodes, 5373 edges.**
Report: **Import Cycles: none.** Probe 1 bounce is gone: EXTRACTED
`shell_graph()` in `blade/shell_graph.py` (L326) does not call
`SectionCut.mid_surface`. `CutWallMetaSDF` is EXTRACTED at
`blade/topology/fields.py` (L98). Canonical loft is `blade/loft/`.

Do not treat these notes as bugs unless a product fails. Do not reopen refused
trailing-edge extract
([SHELL_FE_REPLACEMENT_GAPS.md](../../section/reports/SHELL_FE_REPLACEMENT_GAPS.md)).
Related: [STRUCTURE_PROBES.md](STRUCTURE_PROBES.md),
[PERFORMANCE_WORK_SUMMARY.md](../../section/reports/PERFORMANCE_WORK_SUMMARY.md)
kernel scan.

## EXTRACTED graph

God nodes (`SDFNode`, `BladeGeometry`, `emit_constructive_polygons`,
`MultiCellSection`) are intended hubs. High-betweenness EXTRACTED bridges
(`_emit_section`, `_solid_runs`, `emit_constructive_polygons`) are the station
pipeline. Lazy `cut.py` ↔ `shell_graph.py` wrappers stay; they are not an
EXTRACTED import cycle.

## God nodes

Graphify “god node” is highest **degree** (37 down to 24 edges on this corpus).
It is not a God Class. The top ten are one kernel: CSG base, GEOS
(`MultiPolygon`), loft façade, section host, constructive emit. Do not split
them to lower degree. A problem would be one type that owns unrelated jobs or
cannot be tested without the whole graph.

## `ConstructivePolygons`

`emit_constructive_polygons` returns `dict[str, Polygon | MultiPolygon]`. Keep
the dict. Do not add a dataclass: wrapping `{by_label}` does not move edges off
the builder. Query helpers would duplicate `boundary/labels.py`. C++ analogue is
`std::map<std::string, Polygon>`. Input `section: Any` is a separate DIP item;
do not reopen it here. Do not split walls / cores / cavity unless a caller
needs cavity as a published product (G9 cores stay opt-in).

## Ignore (this filter)

| Report item | Why |
|---|---|
| Surprising connections | **INFERRED** |
| Suggested INFERRED-edge questions | INFERRED |
| Isolated `SectionCut` / `ShellMidlineStrip` | Methods on other node ids |
| Corpus check “4 files” | Incremental detect subset; graph is 3064 nodes |
| `shell_graph()` EXTRACTED-calls `shell_strips_from_graph` at L270 | Same-file AST inversion |

Probes 2 and 5 stay keep. G-TE / G-NET extract stays refused.

## Gaps in this review

Scan is `src/implicit_geometry_kernel` only (not tests, examples, figures).
`--update` left the report corpus stub at four files; community names were
reused. Sixteen percent of edges are INFERRED and were ignored. Lazy imports
often are not EXTRACTED `imports`. The graph is undirected. No pytest, figure
regen, or profile timing in this pass. A full (non-incremental) graphify rebuild
cleans the stub only.

## Product holes this graph cannot close

Refused consumer extract: G-TE … G-ACC. G4 step-`T` locked. WP6 TE classifier
stays an OML-ring helper. G9 solids off the default loop. G10/G11 closed here.
No **Open** BG/SG/CF rows.

## Larger refactors

**Not necessary.** [REFACTOR_PLAN.md](../../REFACTOR_PLAN.md) steps 1–8 done;
step 9 declined. Merged-shell waves 0–7 shipped. STRUCTURE_PROBES 1/3/4 landed.
Do not flatten mixins, merge `StationNetwork`/`ShellGraph`, wrap emit in a class,
or port G-TE/G-NET. C++ maps the same types.

Outstanding in `src/`: none. Next kernel work is a failing figure or test, or a
consumer pin via `scripts/check_consumer.py`.

## Release (not v1.0.0)

Package version stays **0.9.0**. Do **not** tag `v1.0.0`. SemVer 1.0 would freeze
the public API; no consumer has pinned this kernel and `check_consumer.py` is
unused. Programme-complete is not 1.0. A later maintainer pin is **0.9.1**
(docs/alias) or **0.10.0** (minor Unreleased), after that work is committed.
Parked WP6 / locked G4 can still force a break; that is cheaper on 0.9.x.

## Tests

No extra suite to close this review. CI: `pytest -q` (3.11–3.13), import
`core-install`, `scripts/run_examples.py`, ruff, mypy. Golden snapshots in
`tests/test_snapshots.py`. Downstream `check_consumer.py` is unused until a
consumer pins a tag. Do not add a line-% gate or landmark tests that freeze a
consumer extract.
