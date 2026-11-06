# structure

`implicit_geometry_kernel.structure` describes what a thin-walled structure is made of
and combines the parts into signed-distance fields. It is geometry only: thicknesses are
span laws in metres, and plies, materials and meshes belong to downstream packages.

It sits between `sdf` and `section` and imports only `sdf` (enforced by
`tests/test_import_boundary.py`).

## The four contribution kinds

| Kind | What it is | Examples |
|---|---|---|
| `WallContribution` | Part of the wall. One base wall covers the whole OML; a partition sets its own full thickness inside a `footprint`. | skin, spar caps, local doublers |
| `BranchContribution` | An internal member with its own mid-surface that meets the wall. Optional `clip`, standing-off `depth` with a `free_end`, and a span `extent`. | shear webs, stringers |
| `RegionContribution` | A volume with no mid-surface. | sandwich cores, bonds |
| `TransverseContribution` | A member in the section plane at one station `z`, with a spanwise thickness. | ribs, bulkheads, frames |

## Fields of an assembled `Structure`

With `T(x, y, z)` the wall thickness (base law, raised to a partition's law inside its
footprint):

- `wall_mid = phi_OML + T / 2`: zero on the wall reference surface;
- `cavity = phi_OML + T`: negative in the hollow inside the wall;
- `material`: negative inside any member, `meta` = the owning member's feature id. Where
  members overlap the higher `FeatureSpec.priority` wins (core 0, skin 1, web and
  stringer 2, cap 3, transverse 4, bond 5 by default). It is clipped to the OML.

`events` are the stations where a branch starts or ends inside the span;
`forced_stations()` adds the transverse members. A spanwise mesh should contain all of
them so its topology changes only at a station.

## Adding a member type

Subclass `MemberBuilder`, return contributions from `build(ctx, prior)`, and register
it. Builders run in order and see earlier output, so a member that joins others comes
after them. No other code changes. `tests/structure/test_assembler.py` builds a box
beam from all four kinds and adds a new kind this way.

## Span laws

`ConstantLaw`, `TableLaw` (monotone cubic PCHIP or linear, held constant past the end
knots), `ScaledLaw`, `SumLaw`, `ProductLaw`. They are picklable objects, not lambdas.
