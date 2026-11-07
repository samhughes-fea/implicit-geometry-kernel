# structure

`implicit_geometry_kernel.structure` describes what a thin-walled structure is made of
and combines the parts into signed-distance fields. It is geometry only: thicknesses are
span laws in metres. Constitutive plies and meshers belong to consuming tools;
this package supplies the geometry seams (`ThicknessLaw` in metres, named edges).

It sits between `sdf` and `section` and imports only `sdf` (enforced by
`tests/test_import_boundary.py`).

## The four contribution kinds

| Kind | What it is | Examples |
|---|---|---|
| `WallContribution` | Part of the wall. One base wall covers the whole OML. A partition sets the full wall thickness inside a `footprint` (one laminate); a `stacked` layer lies beneath the wall above it. | skins, spar caps, local doublers |
| `BranchContribution` | An internal member with its own mid-surface that meets the wall. Optional `clip`, standing-off `depth` with a `free_end`, and a span `extent`. | shear webs, stringers |
| `RegionContribution` | A volume with no mid-surface. | sandwich cores, bonds |
| `TransverseContribution` | A member in the section plane at one station `z`, with a spanwise thickness. | ribs, bulkheads, frames |

## Fields of an assembled `Structure`

With `T(x, y, z)` the wall thickness (the base law, replaced by the highest-priority
partition inside its footprint, plus any stacked layers):

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

## How sections use it

`MultiCellSection` is assembled through this package (chord frame, before twist):

| Section member | Contribution |
|---|---|
| suction-side skin | base wall |
| pressure-side skin | wall partition on its side of the camber line |
| spar caps | stacked wall layers beneath the skin |
| shear webs | branches whose mid-surface is the web axis |
| sandwich cores, edge inserts | regions |

The builders live in `section/assembly/members.py`; `plan.py` resolves the inputs
first. Each contribution carries the section component it made as its `source`, and the
section's labelled components are those, in build order. `MultiCellSection.structure`
returns the assembled structure; pass `registry=` to add builders.

Caps are a separate layer under the skin, so the geometry matches the section's
existing walls. A merged cap-and-skin laminate (one mid-surface) is planned as an
opt-in section option.
