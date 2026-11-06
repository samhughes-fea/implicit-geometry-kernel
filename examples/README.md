# Examples

Worked examples of the kernel, one folder each with its script, the figures it writes
and a README that walks through the model and the results. Each one checks the kernel
against something independent: a closed form, Shapely, a graph identity or a fresh
integration.

| # | Example | What it shows | Main API |
|---|---|---|---|
| 01 | [Loft and cut](01_loft_and_cut/) | A tapered, twisted wing lofted from two stations and cut anywhere | `BladeGeometry`, `section_at` |
| 02 | [Layout catalogue](02_layout_catalogue/) | All 40 layouts on one twisted airfoil; webs, caps and cells checked against each spec, and the cells against the offset outer mould line | `SYSTEM_TYPE_KEYS`, `resolve_system_type` |
| 03 | [Signed distances and CSG](03_signed_distance_csg/) | A box beam from primitives; which operators keep a distance; grid against quadtree integration of a thin tube; Steiner's formula | `sdf`, `SDFGrid`, `integration` |
| 04 | [Section properties](04_section_properties/) | Wall areas and second moments of the solids against a circular tube's closed forms, and what the midline idealisation costs | `SectionGeometry(outputs=("properties",))` |
| 05 | [Mid-surface and cells](05_mid_surface_and_cells/) | The hand-off to meshing: walls, cells and the merged mid-surface, with Euler's formula, normals and stack offsets checked | `SectionCut.walls`, `cells`, `mid_surface()` |
| 06 | [Web placement](06_web_placement/) | Webs as chord fractions or fixed offsets along a tapered blade; where a straight rear web is dropped, against the airfoil thickness | `web_positions`, `web_offsets`, `dropped_webs` |
| 07 | [Airfoils and root transition](07_airfoils_and_root_transition/) | NACA, ellipse and flatback profiles; a circular root lofted into an airfoil | `Airfoil`, `Station` |
| 08 | [Frames and twist](08_frames_and_twist/) | NACA, S and B frames; second moments rotated with the outline; the sign of twist | `frames` |
| 09 | [Reference blades](09_reference_blades/) | Five windIO 2.0 blades: 2-D cuts, planform, elevation/plan/iso, landmark close-ups; optional Gmsh probe | `BladeGeometry`, `Airfoil` |
| 10 | [Transverse rib](10_transverse_rib/) | A midspan rib as closed faces on the shell graph at that z only | `TransverseMember`, `ShellGraph.faces` |
| 11 | [Partial-span stringer](11_partial_span_stringer/) | A free-end stringer on the graph between two span stations | `Stringer`, `ShellGraph.branches` |

## Running

```bash
pip install -e ".[dev]"
python examples/02_layout_catalogue/layout_catalogue.py   # one example
python scripts/run_examples.py                            # every example and demo, as CI does
```

Each script takes a few seconds and writes its figures into its own `figures/` folder.

## Layout

```
examples/
├── README.md             this index
├── _shared/plotting.py   plot style, wall colours (skin blue, caps red, webs green), draw_cut(), save()
├── NN_<name>/            README.md, <name>.py, figures/*.png
├── section/              short print-only scripts for the section subpackage
└── demos/                figure scripts carried over from blade-structure
```
