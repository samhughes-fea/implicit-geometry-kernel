# R · Reference theory

The R manual derives the geometry the kernel computes. Every formula the code
evaluates is derived here, equation by equation, from a stated definition or a
cited result; the other manuals cite these documents by key ([R1]) and their
equations by number, for example (R1.4), rather than restating them.

## Documents

| Key | Title | Status |
|---|---|---|
| R0 | [Notation and bibliography](R0_notation.md) | Written |
| R1 | [Coordinate frames and twist](R1_frames.md) | Written |
| R2 | Aerofoil parametrisation and families | Planned |
| R3 | Planform and span laws | Planned |
| R4 | Signed distance fields and constructive solid geometry | Current text: [sdf](../sdf/README.md) |
| R5 | Section layout: skin, caps, webs and the WX-Y-Z-NR catalogue | Current text: [section](../section/README.md), [taxonomy](../blade/concepts/taxonomy.md) |
| R6 | Mid-surfaces, cells and the thin-wall idealisation | Current text: [section](../section/README.md), [structure](../structure/README.md) |
| R7 | Section properties by quadtree and shoelace integration | Current text: [integration](../integration/README.md), [shoelace](../integration/shoelace.md) |
| R8 | Spanwise members and the loft | Current text: [blade](../blade/README.md), [structure](../structure/README.md) |
| R9 | Shell mesh generation | Current text: [geometry to mesh](../section/guides/geometry_to_mesh.md) |
| R10 | Feasibility, drops and failure modes | Current text: [feasibility and drops](../section/guides/feasibility_and_drops.md) |

Documents marked "current text" keep their present page until they are rewritten
to the template below; zero-level extraction ([extract](../extract/README.md)) is
folded into R4 and R9 when they are.

## Document template

Each R document has the same eight parts, in this order:

1. **Problem.** One paragraph: what is computed, from what.
2. **Definitions and assumptions.** A numbered list; each assumption states when it
   holds.
3. **Derivation.** One step per numbered equation, from a definition or a cited
   result to the implemented formula. The sentence between two equations names the
   operation that takes one to the next. Every approximation is marked where it
   enters, with its order and range of validity.
4. **Algorithm.** Numbered steps where the computation is more than a formula.
5. **Implementation.** The function that evaluates each result, by name.
6. **Limits.** What the formulas do not cover, and the open defects that affect
   them (bug log IDs).
7. **Verification.** The test or V document that checks each result against a
   closed form or an independent computation.
8. **References.**

```{toctree}
:maxdepth: 1
:caption: R documents

R0_notation
R1_frames
/sdf/README
/section/README
/blade/concepts/taxonomy
/structure/README
/integration/README
/integration/shoelace
/blade/README
/extract/README
```
