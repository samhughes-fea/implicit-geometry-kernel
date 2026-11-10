# Changelog

## Unreleased

- CSG property panels reuse occupancy thin-band wall-tube refine (BG-025).
  Occupancy ``max_cells`` is 400000 so 5D root stations finish (BG-026).

- Docs hygiene: package bug logs (no repo-root ``BL-`` file); section module
  table uses ``walls/`` / ``boundary/`` / ``properties/``; reports index;
  span interpolation is linear (PCHIP wording dropped from the refactor plan).

- House-rule splits: ``blade/loft/`` holds the loft façade; ``windio`` nd-arc
  helpers are ``windio_nd.py``. ``members.py``, ``mesh/polygon``, ``spec/runner``,
  ``quadtree/soa``, and ``plot_section_hybrid`` stay under 400 lines. Public
  names are unchanged. ``build_spar_caps`` / ``skin_split`` remain under 400;
  station network and spanwise SDF fields live in ``blade/topology/``; planform
  CSV I/O in ``blade/concepts/planform.py``.

- Teaching scripts, section demos, contract 01, feature-field loft and box-beam
  assemble use ``SectionSpec`` / ``SectionGeometry(spec, options)`` and
  ``TagNode(Circle|Polygon, …)``.

- Manuals show ``SectionGeometry(spec, options)``; glossary topology homes are
  ``blade/topology/``. ``graphify-out/`` and local GEOS scratch are gitignored.

- ``SectionGeometry(spec, options)`` is the constructor. Keyword
  ``layout=`` / ``airfoil=`` and ``from_kwargs`` warn and are removed in the
  next minor. ``BladeGeometry.section_at`` builds ``SectionSpec`` /
  ``SectionOptions``; min-thickness scale-up is ``_min_thickness_scale``.
  ``BladeGeometry.layout`` is typed ``str | SystemLayoutSpec | SectionLayout``.

- ``TagNode`` is the only meta stamper. ``ExprSDF``, ``CircleSDF`` and
  ``PolygonSDF`` wrap ``TagNode`` and warn; removed in the next minor.

- Spar-cap construction dispatches ``CapPlacementKind`` through a callable
  table in ``build_spar_caps``.

- ``WallRegionGraph`` stores ``cap_placement`` (``CapPlacementKind``).
  ``cap_mode`` (``band``/``pads``/``box``) is a deprecated alias. Topology
  compile lives in ``blade/topology/layout_compile.py``; loop signatures in
  ``signatures.py``.

- Rib mesh patches split only at live web feet (``n_webs + 1``). Cap
  ``LoopBranch`` breakpoints stay on the skin midsurface and no longer seed
  extra rings (Family B/C).

- Example 09 lofts 18 catalogue tokens on the scaled 12 m planform (adds
  ``0A-LC-0R``, ``0B-LC-0R``, ``1B-F-LC-10R``, ``2A-F-LC-10R``, ``2C-F-LC-10R``,
  ``2D-CN-LC-10R``).

- Example 09 synthetic planform is the original 20 m knots scaled to 12 m
  (span, chord, and wall gauges × 12/20), not a truncated tip.

- Meshing lives under ``mesh/``: kernel polygonal fill (``mesh/polygon.py``,
  GEOS Voronoi), numpy Steiner Delaunay (``mesh/delaunay.py``), and Gmsh in
  ``mesh/gmsh/`` (  ``plane.py`` fills all rib cells in one Gmsh session). Native rib fill uses
  scipy Qhull when present, else numba Bowyer–Watson. Rib n-gon fill
  (``mesh/polygon.py``) uses Qhull Voronoi plus a numba clip; faces have more
  than four nodes. Example 08/09 rib figures use that n-gon fill. OBJ still
  writes the cell n-gon.

- Rib drawing/STL ear-clips each cell n-gon (same loft nodes, no Steiner
  points) so matplotlib no longer fans through the LE cell. Rib mesh strokes
  match loft quads on webs and skins; cutaway captions note the omitted-suction
  camera.

- Example loft palette ``RETRO_PANELS`` is 32 colours; 08/09 loft quads share
  a darker ``MESH_EDGE`` stroke.

- Breaking: PA web freeze uses a uniform root-band lattice (up to three
  interiors of ``[x_LE, x_TE]``); extra webs continue that lattice (first extra
  at root TE), replacing the authored ``_PA_SLOTS`` fractions. Example 09 adds
  ``5D-F-PA-10R`` (second secondary at two ``Δx``).

- Breaking: span interpolation is linear only. ``TableLaw(kind=)`` and
  ``BladeGeometry(span_kind=)`` are removed. Ply-drop gauges are plateaus with
  linear ramps (``TableLaw`` / ``PlyDrop``); cubic PCHIP is gone.

- ``SectionCut.layout`` is the authored blade concept after web drop, not a
  sibling catalogue token. Realised W and X are ``station_webs`` /
  ``station_family``. Cascade helpers remain assemble recipes for that slice.

- Docs and gallery cleanup after the closed architecture programme: stale
  Shapely / Phase 2d / CAD-wait wording in maintainer notes; ``[mesh]`` dropped
  from the public package map; glossary rows for ``station_webs`` /
  ``station_family``; ``examples/demos/section`` hosts assemble through
  ``SectionGeometry``.

## 0.8.2

- Family **C** pad 0 is stored as ``spar_cap_*_0`` in the unrotated pad dict.
  Constructive union is indexed pads only; the canonical label stays the merged
  product (and family **D** box).
- Mapped loft strips resample unmatched station polylines to a shared vertex
  count and omit collapsed or non-finite faces (example 09 3-D autoscale).

## 0.8.1

- Family **C** constructive cap union includes pad 0 (stored under the
  canonical surface label). Hybrid ``2C-CN-LC-0R`` ``spar_cap_suction_surface``
  area 9.15e-4 → 1.83e-3 (×2.00); pressure 9.13e-4 → 1.82e-3. Cores lose the
  extra pad. Quadtree occupancy was already the pad union; hybrid parity for
  C-family wall solids uses the same tight tolerance as the other families.

## 0.8.0

- Breaking: ``apply_shapely`` is removed. Use ``apply_clip``. The name was a
  deprecated alias in 0.7.0.

## 0.7.0

- CI after the tree freeze: mypy on loft mixins and GEOS ctypes; named
  ``feature_curve`` passes all four ``partition_edge_curve`` arguments;
  ``PolygonOps.distance`` is GEOS ``GEOSDistance_r`` (tests use ``OPS``, not
  methods on kernel ``Polygon``); C-family hybrid vs quadtree area is
  documented as per-web pads, not circular disks; ``09_blade_concepts``
  ignores empty/non-finite 3-D faces and degenerate axis limits. Public API
  snapshot records the GEOS / ``apply_clip`` names; ``apply_shapely`` is a
  deprecated alias of ``apply_clip``.

- Shapely is gone. Constructive offset and overlay stay **GEOS** (JTS C++)
  through ``implicit_geometry_kernel.geos`` (ctypes ``libgeos_c``, kernel
  ``Polygon``, WKB). Install system/conda GEOS or set ``GEOS_LIBRARY_PATH``.
  ``GeosPolygonOps`` is production ``OPS``. Emitters and ``apply_clip`` route
  Booleans through ``OPS``. Catalogue section goldens follow the kernel WKB
  ``Polygon`` (same GEOS as Shapely 2 on this machine; previous JSON was
  slightly off).

- Mid-surface patches: STL/OBJ per authored member; STEP via extra ``[cad]``
  (cadquery compound, no fuse). ``members=`` filters patches. Debug SDF contours
  declare extra ``[contours]`` (contourpy).

- Teaching example ``06_consumer_seams.py`` shows how a laminate or
  section-model repo consumes ``ThicknessLaw``, ``named_edges``, ``size_fields``
  and ``merge_cap_and_skin``.

- Architecture programme closed at the target tree: ``catalogue/`` is a layer;
  ``plotting``, ``mesh``, ``blade.cad``, and ``profiling`` are extras. House-rule
  splits: loft façade, topology map, shell-graph build, spar-cap pads, section
  plotters. ``StiffenerWall`` is a registered constructive kind; ``build`` is the
  only ``WallType`` pipeline hook. Extract tests live under ``tests/extract/``.

- Constructive polygon calls go through ``PolygonOps`` (gateway
  ``section/boundary/polygons.py``). ``apply_clip`` is the clip interpreter;
  ``apply_shapely`` is a deprecated alias (removed in 0.8.0). Shoelace and
  plotting do not import Shapely. Oracles live in
  ``tests/section/test_polygon_ops_oracle.py``.

- Target tree lists ``catalogue/`` as a layer; ``plotting``, ``mesh``, ``blade.cad``,
  and ``profiling`` are extras. Import-boundary tests cover ``mesh``. ``SDFGrid``
  debug contours use contourpy only (no matplotlib fallback).

- Station wall kinds dispatch on ``SECTION_BUILDERS`` flags (``on_shell_product``,
  ``polygon_by_surface``). ``build_section_view`` is the assemble host helper in
  ``section.assembly``; gated products stay on ``build_section_geometry``.
  Catalogue and concept shims (``section.layout.system_layout``, ``blade.concept``)
  are gone. ``blade.topology`` is a package; OML writers live in ``blade.oml``.

- ``sdf`` is one tree: ``kernels``, ``nodes``, ``shapes``, ``solve``. Retired
  top-level ``sdf/{ir,airfoil,csg,grid,primitives,naca,transforms,_polygon_numba}``.
  One ``AirfoilSDF`` (``sdf.shapes``); author via ``Airfoil.to_sdf()``.
  ``SDFGrid`` is an evaluation grid; product integrals stay in ``integration``.
  ``GeometrySolver`` stays in ``sdf.solve``; ``extract`` traces zero-set polylines.

- Catalogue tokens live in ``implicit_geometry_kernel.catalogue`` (below
  ``section``). ``blade.concepts`` authors a full blade (``BladeConcept.from_key``,
  R-letter ribs). Taxonomy: ``docs/blade/concepts/taxonomy.md``.

- Documentation is split into a **manual** (package READMEs, guides, frames,
  taxonomy, numbered equations) and **internal** pages (bug logs, ADRs,
  reports). House style: ``docs/conventions/DOCUMENTATION_STANDARDS.md``.
  Completeness: ``docs/glossary.md`` lists every kernel ``__all__`` name.
  Hub data flow is ``BladeGeometry`` → ``SectionCut`` → ``SectionGeometry``.
  ``docs/contract/README.md`` describes live DTOs only (loads remain in
  ``docs/FORK.md``). Sphinx/autodoc is deferred.

- Mapped loft and geometry-to-mesh handoff are frozen (consume table in
  ``docs/section/guides/geometry_to_mesh.md``; loft rules and closed bugs in
  ``docs/blade/bug_log.md``).

- Loft keep-set equals intercept ``DroppedWeb.index`` when that set is passed:
  no 0.08 cluster or extra TE-plug omit. ``_regularise_anchors`` keeps a
  ``web_*`` over a nearby ``cap_*``. TE wrap walks include the OML ridge;
  TE metrics use nearest-ridge / max chord-frac, not blade-x.

- TE wrap is a loft-zone event from the skin-mid opening (blunt vs sharp).
  Blunt stations walk the real TE edge; sharp stations share one ridge key
  and drop wrap. Hysteresis holds the boolean so a station at ε cannot
  chatter. `_te_frac_drop` is superseded (BG-005).

- Loft VTK writes PNG ``panel_colour`` RGB as ``COLOR_SCALARS`` so ParaView
  shows the same patches; cell field ``role`` keeps the region name.

- Skin loft at a PA web onset/drop rules each t1 panel from the t0 mould walk
  between its mapped ends (not a cavity chord or a collapsed TE point). Dropped
  web anchors leave the loop.

- PA/LE loft anchors follow the metre freeze, not catalogue LC fractions, so
  ``web_0`` / ``web_3`` walk to the LE/TE with the live intercept instead of
  vanishing mid-chord.

- Cylinder loft clustering no longer uses a 0.20 chord radius, so a 4-web PA
  lattice (0.2 / 0.4 / 0.6 / 0.8) keeps ``web_2`` through the root. TE-named
  plugs are still omitted on a cylinder.

- Loft templates keep authored web indices from ``DroppedWeb.index`` (order-
  preserving match only as fallback). TE corners track the OML TE so cylinder–
  foil stations do not jump. Example 09 tip chord is 0.75 m so 3D-F still
  anchors one web at the tip.

- Example loft dumps write OBJ, VTK POLYDATA, and NPZ (`polygons` includes rib
  n-gons). Example 09 also writes a suction-omitted `cutaway.obj` / `cutaway.vtk`.

- Ribs (and tip closure) are n-gons whose vertices are the mapped-loft skin
  and web nodes at the member station. Skins and webs stay span quads. OBJ
  keeps those n-gons; STL fans from the first loft vertex. TE keys are the
  start of each TE→LE side arc.

- Breaking: PA web freeze is hard-coded blade-x from the pitch axis, not a
  ``k·c_root/n`` lattice and not the two-web ``(0.15, 0.50)`` pair. Primary
  slots sit in the root ``[x_LE, x_TE]`` band about PA; secondary slots sit in
  the TE shelf ``[x_TE(root), x_TE(max_chord)]``. LC/LE remain authored chord
  fractions / metres from LE. ``BladeConcept`` lives in
  ``blade.concepts`` (``blade.concept`` is a shim). Web labels keep authored
  indices through drop/onset. Example 09 adds ``4D-F-PA-10R``.

- Breaking: `BladeGeometry` requires `pivot_fraction` (float or span table) —
  the S-frame origin and XY stacking point. There is no 0.25 / 0.333 fallback.
  `pitch_axis=` remains an alias for the same table; do not pass both.

- Authored loft topology: `WallRegionGraph` (concept `from_section_layout` and
  windIO `from_windio_structure`), `ShellTopologyTemplate.map_vertex`,
  `zone_region_graphs`, and `loft_mapped_strips`. Figure loft and
  `write_mid_surfaces` join by region id, not run-name intersection. windIO
  compile adds a `tip_closure` `TransverseMember` at last span (plus YAML
  `ribs` / `bulkheads`). Example 08 colours loft faces with `RETRO_PANELS` by
  role and draws rib faces. Cores stay solids. `shared_loft_names` is unchanged
  (Wave 2). No Gmsh / PCOMP.

- Locked loft: `shell_loft_zones`, `LockedLoftPlan`, and `resample_locked` keep
  held chord counts and shared partition vertices. Surviving names loft through
  topology events; only missing names are omitted. `write_mid_surfaces` and
  example 08 use that rule instead of a whole-blade name intersection.

- Conformal loop: cap paint on two web-bounded panels stays two runs
  (no jump-join). Published `chord_window` survives a web drop instead of
  renaming caps to `per_web_pads` `*_0`.

- windIO `cap_windows`: `CapPlacementRule.from_chord_windows` stores the
  published `(lo, hi)` chord interval (`chord_window`). Assemble uses that
  band, not a pitch-centred half-width pad.

- Shell mesh: conformal **skin mid** (no stack-NA kink at cap ends). Cap
  footprints are `LaminateOffsetRegion`s (`ShellGraph.offset_regions()`) for
  downstream constitutive OFFSET. `wall_mid_sdf` follows the skin mid; stack
  `T` stays a thickness oracle.

- Breaking: `result.midline_strips` / `build_shell_midline_strips` are
  skin+web only. Cap solids stay on `cut.walls` (and properties); they are
  not a second mesh polyline. Anyone counting cap rows on the midlines
  product must use `cut.walls` or `shell_graph` partitions.

- windIO: span-varying `web_position_table_from_windio_structure` /
  `cap_window_table_from_windio_structure`. Example 08 compiles published
  `nd_arc` tables instead of overlaying a 2D-F layout.

- Docs: layered `README.md` files from the repository root through `docs/`,
  `src/`, `tests/`, and `examples/`. New hub at `docs/README.md` with a
  capability map and mermaid package / data-flow trees. Package pages for
  `blade` and `sdf`; restored `docs/conventions/` (documentation, API runner,
  module layers).

- Breaking for scripts and bookmarks: numbered examples now live under
  `examples/<package>/NN_*` (blade, section, structure, sdf, frames). Teaching
  scripts moved to `examples/section/teaching/` including
  `05_midline_handoff.py` (replaces the deleted mesh-driver twin).

- WallType: leftover label tests go through `wall_type_for_label` / `kind`
  (cut key-points, zone gauges, web counts, property wall-solid keys). A
  `BondlineBuilder` is registered on `SECTION_BUILDERS` with an empty `build`
  so catalogue layouts stay at 40 letters; tests may override it.

- G9: `shell_graph(..., include_solids=True)` adds core and insert `ShellRun`s
  from polygon boundaries (`surface` `"core"` / `"insert"`). Default remains
  off; 2D-F run-name goldens unchanged.

- WP6: named `feature_curve(geo, name_a, name_b, zs)` delegates to
  `partition_edge_curve`. `classify_trailing_edge` accepts a stacked skin loop
  (`ShellGraph`) as well as an OML ring. `LoftedPolygonSDF` caches an AABB /
  uniform-grid BVH for batched queries; `phi` is unchanged.

- `BladeGeometry.write_mid_surfaces` writes one STL/OBJ patch per loftable
  run. Optional extra `[cad]` (cadquery) exports those patches as STEP from
  `blade.cad` only.

- Additive `SpanwiseScheme` and `BladeConcept` compile to existing
  `BladeGeometry` kwargs (`from_concept` / `concept`). Catalogue keys stay
  presets of `SectionLayout`.

- `sdf` remains `implicit_geometry_kernel.sdf`. Import tests enforce
  down-only layers (`contract` ← `sdf` / `integration` / `frames` ←
  `structure` ← `section` ← `blade` ← `extract`).

## 0.6.0

- Breaking: 0.5.x deprecation aliases are gone. Removed dual names
  (`CellLoop`, `SpanWeb`, `SectionGeometryRunner`, `station_at`,
  `web_x_positions_m`), shim packages (`section.api`, `section.placement`,
  `section.constructive`, `section.props`, `section.recipes`, root `plotting`,
  flat `sdf.*` redirects, `integration.quadrature`), the keyword form of
  `build_section_geometry`, `compile_expr_to_node` / `eval_via_nodes`, silent
  `*Node` CSG aliases, and `blade.cut` re-exports of `shell_graph`. Import
  `implicit_geometry_kernel.plotting`, `section.spec`, and
  `blade.shell_graph`. Geometry goldens unchanged.

- Architecture: this kernel is the source of truth for blade geometry
  parameterisation. Plots go behind `[plots]`; mid-surface field from
  `shell_graph`; demos stay and will be ported; seams for thickness/edges/size
  live here; no Gmsh extra, no constitutive plies, no separate blade-model
  package. Optional `[cad]` waits until after those seams.

- Plots live in `implicit_geometry_kernel.plotting` behind extra `[plots]`.
  Geometry modules still do not import matplotlib.

- Exactness: `wall_mid_sdf` is the signed polygon SDF of the `shell_graph`
  skin loop at each `z` (run interiors and partition vertices on `phi = 0`).
  Webs remain branches. Section and blade goldens unchanged. Worst `|phi|` at
  loop-run interiors on all 40 layouts is below 5e-4.

- Demos under `examples/demos/` import `section.spec` and
  `implicit_geometry_kernel.plotting`. Anatomy demo also prints a
  `BladeGeometry.shell_graph()` inventory.

- Seams for meshing and thickness: `ThicknessLaw` is `T(z, s)` in metres
  (`from_span` wraps a span law). `BladeGeometry` exposes `thickness_law()`,
  `named_edges(z)`, `size_fields(points)`, ply-drop ends on `forced_stations()`,
  and opt-in `merge_cap_and_skin` (default off; goldens use stacked caps).

- Blade completeness: ply-drop `T(z, s)` ramps, named key points, watertight
  OML STL/OBJ, batched point queries, planform CSV in / windIO dict out, and
  `validate()` with station and member named. NACA 4-digit thickness at the
  nose matches the analytic formula to 1e-8; insert polygons match shoelace
  area to 1e-12. Geometry goldens unchanged.

- Optional `BladeGeometry(span_kind="pchip")` interpolates chord, twist, pitch
  axis and offset with a monotone cubic (airfoils stay linear).
  `Airfoil.cylinder_to_airfoil` is the named root cylinder→airfoil law.

- `Airfoil` grows CST, NACA 5/6, camber-line t/c families, flatback TE,
  Selig/Lednicer readers and interpolation by thickness. NACA vertex maths
  stays in `sdf.shapes.naca`.

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

- Dual names from 0.5.x are removed (see Breaking above).

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
