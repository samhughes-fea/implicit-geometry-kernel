# blade

Airfoils, stations, loft, cuts, shell graph.

**Import:** `implicit_geometry_kernel` (root barrel) or `.blade`.

↑ [package](../README.md) · [docs](../../../../docs/blade/README.md) ·
[examples](../../../../examples/blade/README.md) · [tests](../../../../tests/blade/README.md)

Key types: `Airfoil`, `Station`, `BladeGeometry`, `BladeConcept`, `SectionCut`, `shell_graph`,
`LockedLoftPlan`, `resample_locked`, `shell_loft_zones`.
`section_at(z)` calls `SectionGeometry` and returns a `SectionCut`.
In-memory object graph: [docs/blade](../../../../docs/blade/README.md#in-memory-objects).
Names: [glossary](../../../../docs/glossary.md).
Concepts (`BladeConcept`, LC/LE/PA web positions): [concepts](concepts/).
