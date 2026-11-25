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
| R2 | [Aerofoil parametrisation and families](R2_aerofoils.md) | Written |
| R3 | [Planform and span laws](R3_span_laws.md) | Written |
| R4 | [Signed distance fields and constructive solid geometry](R4_sdf_csg.md) | Written |
| R5 | [Section layout: skin, caps, webs and the WX-Y-Z-NR catalogue](R5_section_layout.md) | Written |
| R6 | [Mid-surfaces, cells and the thin-wall idealisation](R6_midsurfaces_cells.md) | Written |
| R7 | [Section properties by quadtree and shoelace integration](R7_section_properties.md) | Written |
| R8 | [Spanwise members and the loft](R8_loft_members.md) | Written |
| R9 | Shell mesh generation | Current text: [geometry to mesh](../section/guides/geometry_to_mesh.md) |
| R10 | Feasibility, drops and failure modes | Current text: [feasibility and drops](../section/guides/feasibility_and_drops.md) |

Documents marked "current text" keep their present page until they are rewritten
to the template below; zero-level extraction ([extract](../extract/README.md)) is
derived in R4, Section 3.11.

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
R2_aerofoils
R3_span_laws
R4_sdf_csg
R5_section_layout
R6_midsurfaces_cells
R7_section_properties
R8_loft_members
/sdf/README
/section/README
/blade/concepts/taxonomy
/structure/README
/integration/README
/integration/shoelace
/blade/README
/extract/README
```
