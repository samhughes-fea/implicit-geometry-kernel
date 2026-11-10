# Structure

↑ [docs](../README.md) · [src](../../src/implicit_geometry_kernel/structure/README.md) ·
[examples](../../examples/structure/README.md) · [tests](../../tests/structure/README.md)

This package describes the members of a thin-walled structure and combines them
into signed-distance fields. Thicknesses are span laws in metres. Constitutive
plies and meshers belong to consuming tools.

It sits between `sdf` and `section` and imports only `sdf`
(`tests/test_import_boundary.py`).

Let \(\phi_{\mathrm{OML}}\) be the outer-mould signed distance and \(T(x,y,z)\)
the wall thickness (base law, replaced by the highest-priority partition inside
its footprint, plus stacked layers). Then

$$
\phi_{\mathrm{mid}} = \phi_{\mathrm{OML}} + T/2,
\qquad
\phi_{\mathrm{cavity}} = \phi_{\mathrm{OML}} + T.
\tag{1}
$$

The mid-surface \(\phi_{\mathrm{mid}}=0\) is the wall reference. The cavity is
negative in the hollow interior. The material field is negative inside any
member; overlapping members resolve by `FeatureSpec.priority`.

## The four contribution kinds

| Kind | What it is | Examples |
|---|---|---|
| `WallContribution` | Part of the wall. One base wall covers the whole OML. A partition sets the full wall thickness inside a `footprint` (one laminate); a `stacked` layer lies beneath the wall above it. | skins, spar caps, local doublers |
| `BranchContribution` | An internal member with its own mid-surface that meets the wall. Optional `clip`, standing-off `depth` with a `free_end`, and a span `extent`. | shear webs, stringers |
| `RegionContribution` | A volume with no mid-surface. | sandwich cores, bonds |
| `TransverseContribution` | A member in the section plane at one station `z`, with a spanwise thickness. | ribs, bulkheads, frames |

## Fields of an assembled `Structure`

- `wall_mid`: \(\phi_{\mathrm{mid}}\) in \((1)\);
- `cavity`: \(\phi_{\mathrm{cavity}}\) in \((1)\);
- `material`: negative inside any member; `meta` is the owning feature id.
  Default priorities: core \(0\), skin \(1\), web and stringer \(2\), cap \(3\),
  transverse \(4\), bond \(5\). The field is clipped to the OML.

`events` are stations where a branch starts or ends; `forced_stations()` adds
transverse members. A spanwise mesh should contain all of them so that topology
changes only at a station.

## Adding a member type

Subclass `MemberBuilder`, return contributions from `build(ctx, prior)`, and register
it. Builders run in order and see earlier output, so a member that joins others comes
after them. No other code changes. `tests/structure/test_assembler.py` builds a box
beam from all four kinds and adds a new kind this way.

## Span laws

`ConstantLaw`, `TableLaw` (linear interpolation, held constant past the end knots),
`ScaledLaw`, `SumLaw`, `ProductLaw`. They are picklable objects, not lambdas.

## How sections use it

`MultiCellSection` is assembled through this package (chord frame, before twist):

| Section member | Contribution |
|---|---|
| suction-side skin | base wall |
| pressure-side skin | wall partition on its side of the camber line |
| spar caps | stacked wall layers beneath the skin |
| shear webs | branches whose mid-surface is the web axis |
| sandwich cores, edge inserts | regions |

Two hosts remain, and they are not merged. `structure.MemberBuilder.build` is the
SDF host. `SECTION_BUILDERS` in `section/assembly/members.py` is the constructive
host: each builder is a `MemberBuilder` that also carries `WallTypeMixin` flags so
later stages (boundary, midlines, cells, properties) dispatch without name prefixes.
`plan.py` resolves inputs first. Each contribution carries the section component it
made as its `source`. `MultiCellSection.structure` returns the assembled structure;
pass `registry=` to add builders.

Caps are a separate stacked layer under the skin by default (`merge_cap_and_skin=False`),
so golden numbers stay put. With the flag on, `BladeGeometry.section_at` replaces the
stacked walls with one laminate per `shell_graph` run (combined thickness).

## Seams for consumers (geometry only)

These live on `BladeGeometry` / `shell_graph`. Constitutive plies and Gmsh stay out.

| Seam | What the kernel hands off |
|------|---------------------------|
| Thickness | `ThicknessLaw`: \(T(z,s)\) in metres (`structure.laws`; `BladeGeometry.thickness_law`) |
| Edges | named partition and ply-drop curves (`shell_graph.partition_edges`, `named_edges`) plus `forced_stations` |
| Size | `size_fields`: distance, thickness, curvature, member id at sample points |
