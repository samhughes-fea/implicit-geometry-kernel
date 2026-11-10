# Blade structure probes

Dated 2026-11-10. Filled 2026-11-10. Aliases 1, 3, 4 landed the same day.
Investigation plus those aliases. No further restructure in this note.

Five SOLID questions on the loft cut and shell graph. Counts from `src/`,
`tests/`, and `examples/` (not `docs/_build`). Do not treat these probes as
bugs unless a product fails. Do not reopen refused trailing-edge extract or
bisector station networks
([SHELL_FE_REPLACEMENT_GAPS.md](../../section/reports/SHELL_FE_REPLACEMENT_GAPS.md)
source of truth).

Canonical modules are under `blade/loft/` and `blade/topology/`.
Graphify EXTRACTED follow-up (no further kernel step):
[GRAPHIFY_EXTRACTED_NEXT.md](GRAPHIFY_EXTRACTED_NEXT.md).

## 1. SRP — `SectionCut` snapshot vs graph builders

[`shell_graph(cut)`](../../../src/implicit_geometry_kernel/blade/shell_graph.py)
(L326) calls `_mid_surface(cut)` (L333).
[`SectionCut.shell_graph`](../../../src/implicit_geometry_kernel/blade/cut.py)
(L154–161) and `.mid_surface` (L123–137) lazy-import those builders.
[`shell_strips_from_graph`](../../../src/implicit_geometry_kernel/blade/shell_graph.py)
(L270) uses `shell_graph(cut).runs()`.

**Findings.** Import cycle remains: `cut.py` lazy-imports `shell_graph` /
`_mid_surface`. `shell_graph.py` TYPE_CHECKING-imports `SectionCut`. Builders
do not bounce through the methods. `_mid_surface` is defined in
`shell_graph_build.py`.

`cut.shell_graph()` in `src/` (method wrappers, intentional): `oml.py`,
`loft/loft_fields.py` (`BladeGeometry.shell_graph(z)` and related),
`topology/fields.py`, `topology/map_walk.py` (2), `topology/ribs.py` (3),
`extract/curves.py`, `mesh/gmsh/shell.py`, `partition_edge_curve` (L353).
`BladeGeometry.shell_graph(z)` in `loft_fields.py` is a different method
(it calls `section_at(z).shell_graph()`).

`shell_graph(cut)` function: definition, `shell_strips_from_graph`,
`topology/network.py`, `topology/map.py` (`station_key_points`).

`cut.mid_surface()`: tests `test_loft.py` / `test_shell_graph.py`, example
`section/03_mid_surface_and_cells`.

**Recommendation.** Alias. Keep `cut.shell_graph()` / `cut.mid_surface()` as
one-line wrappers. Builders must call `_mid_surface(cut)` / `shell_graph(cut)`,
not bounce through the methods. Lazy imports stay only if the
`cut` ↔ `shell_graph` import cycle remains. C++ analogue: struct methods that
forward to free functions.

**Later.** Alias (done 2026-11-10: builders call `_mid_surface(cut)` / `shell_graph(cut)`; methods remain wrappers).

## 2. ISP — `StationNetwork` vs `ShellGraph`

[`station_network(cut)`](../../../src/implicit_geometry_kernel/blade/topology/network.py)
(L59–70) wraps `shell_graph(cut).runs()` in `WallBranch` plus `FeatureTable`.
The package barrel exports both `ShellGraph` and `StationNetwork`.
`station_network_from_meta` remains a meta split (`thickness=0.0`).

**Findings.** `station_network(cut)`: `topology/network.py`, tests
`test_station_network.py` and `test_shell_graph.py`, example
`11_feature_station_network`. `station_network_from_meta`: same example,
`tests/extract/test_wall_meta_extract.py`, `tests/test_shell_graph_gaps.py`
(D1.3 thickness 0). Root and `blade` barrels export both types. `ShellGraph`
has a much larger test surface (`test_shell_graph.py` and related).

**Recommendation.** Keep both. `StationNetwork` is the FeatureTable-tagged
mesher view. Do not merge types. `station_network_from_meta` stays a meta
split.

**Later.** Keep both.

## 3. DIP — `SectionCut.host`

Field at [`cut.py` L104](../../../src/implicit_geometry_kernel/blade/cut.py):
`host: SectionMidlineHost | None`. Bond `host_label` is a different name.

**Findings.** Loft assigns `host` as `result.section` in
[`loft/loft_cut.py`](../../../src/implicit_geometry_kernel/blade/loft/loft_cut.py)
L197. `SectionGeometryResult.section` is
[`SectionMidlineHost`](../../../src/implicit_geometry_kernel/section/midlines/host.py)
(`section/spec/results.py` L145). Reader:
`examples/_shared/blade_viz.py` (`_sample_property_qt` requires `cut.host`,
then `SectionFieldEvaluator.from_section`). Test:
`test_cut_host_is_section_midline_host` in `tests/test_loft.py`.

**Recommendation.** Type the field as `SectionMidlineHost | None`. Do not drop
it; viz depends on it.

**Later.** Type as `SectionMidlineHost | None` (done 2026-11-10).

## 4. SRP — `key_points` vs `place_station_keys`

[`station_key_points`](../../../src/implicit_geometry_kernel/blade/topology/map.py)
(L41) builds LE/TE, cap-edge vertices and web feet via `shell_graph(cut)` and
`wall_type_for_label`. [`SectionCut.key_points`](../../../src/implicit_geometry_kernel/blade/cut.py)
(L163–167) is a one-line call. [`place_station_keys`](../../../src/implicit_geometry_kernel/blade/topology/map.py)
(L65) takes cut, graph, region graph, and `z`.

**Findings.** `SectionCut.key_points` / `BladeGeometry.key_points(z)` (mixin
in `loft_fields.py` L122): example `06_blade_queries`, test
`test_blade_completeness.py`. `place_station_keys`: `topology/map.py` (also
called from itself for wrap), `extract/curves.py`, tests
`test_shell_topology.py` and `test_sizing.py`. Two APIs: LE/TE plus cap/web
feet on the cut vs named loft keys on the graph.

**Recommendation.** Alias. Keep `SectionCut.key_points` as a thin wrapper.
Move the LE/TE/cap/web-foot logic next to `place_station_keys` in `topology/`.

**Later.** Alias to topology (done 2026-11-10: `station_key_points` in `topology/map.py`).

## 5. SRP / composition — `BladeGeometry` mixins

[`BladeGeometry(LoftCutMixin, LoftFieldsMixin, LoftLayoutMixin)`](../../../src/implicit_geometry_kernel/blade/loft/loft.py)
(L73). After the class body, `loft.py` L292–297 copies public callables from
`LoftCutMixin` and `LoftFieldsMixin` onto `BladeGeometry`. Mixins are not
imported as public types from `blade/__init__.py`.

**Findings.** Three mixin modules exist only to split `BladeGeometry`.
`LoftLayoutMixin` is inherited; cut and fields are also copied with
`setattr`. No other subclass of those mixins in `src/`.

**Recommendation.** Keep. File split of one type (C++ analogue: one class,
methods in several translation units). Do not add extra public types.

**Later.** Keep mixins.

## Out of scope here

A consumer's own extract (bisector, MERGED TE, continuous `T`).
A bug-log row (no product failure in these probes).
Further rewrites of `cut.shell_graph()` call sites that are not builders.
