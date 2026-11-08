# 09 · Blade concepts

Catalogue layouts (`BladeConcept` / `SYSTEM_TYPE_KEYS`) lofted on a **made-up**
planform. This is not a published turbine: airfoils, chord, and twist are invented.
Published windIO blades stay in [08_reference_blades](../08_reference_blades/).

```bash
python examples/blade/09_blade_concepts/blade_concepts.py
```

Eleven tokens share one 20 m PCHIP loft (circular root → NACA family, eleven
design knots, chord 1.8 → 2.4 → 0.45 m, twist 12° → −2°).

## Web anchoring

Catalogue tokens: [WX-Y-Z-NR taxonomy](../../../docs/blade/concepts/taxonomy.md).
Code: `implicit_geometry_kernel.blade.concepts`.

All three modes use **hard-coded positions** in native units. `n_webs` picks how
many authored slots to use; it does not invent a lattice.

| Mode | Quantity | Planform |
|---|---|---|
| `LC` | chord fractions | slide with local chord (not parallel to PA) |
| `LE` | metres from the LE (root chord as the stick) | track the LE / taper (not parallel to PA) |
| `PA` | blade-x from the pitch axis | **parallel to PA** — primary = root `[x_LE, x_TE]` about PA; secondary = `[x_TE(root), x_TE(max chord)]` |

Primary/secondary **planform bands are PA-only**. LE/PA omit a web that misses
the inner skin (the same index may seed later). Rib count is required: `0R`
none, `4R` four faces, `10R` ten faces from just after the root through the tip.

| Key | Layout |
|---|---|
| `2D-F-LC-0R` | two webs, box cap, local chord, no ribs |
| `2D-F-LC-4R` | same, four ribs |
| `2D-F-LC-10R` | same, ten ribs |
| `2D-F-LE-0R` | two webs, box cap, leading-edge offsets, no ribs |
| `2D-F-LE-4R` | same, four ribs |
| `2D-F-LE-10R` | same, ten ribs |
| `2D-F-PA-10R` | two webs, box cap, PA blade-x (primary about PA), ten ribs |
| `3D-F-LC-10R` | three webs, box cap, local chord, ten ribs |
| `3D-F-LE-10R` | three webs, box cap, leading-edge offsets, ten ribs |
| `3D-F-PA-10R` | three webs, box cap, PA blade-x (primary about PA), ten ribs |
| `4D-F-PA-10R` | four webs, box cap, PA blade-x (primary + TE-shelf secondary), ten ribs |

Each concept compiles `BladeConcept(key)` → `BladeGeometry.from_concept`. Wall strips use the authored `WallRegionGraph` and `loft_mapped_strips`. Web ends are the shared skin-loop keys (not T-joint feet). Unnamed 3-D PNGs use the leading-edge camera; `*_te.png` is the trailing-edge twin.

Figures per key under `figures/<key>/`:

```
cross_sections/sections.png
regions/{root,transition,max_chord,span70,tip}.png
regions/{root,transition,max_chord,span70,tip}_te.png
region_cutaways/{root,transition,max_chord,span70,tip}.png
region_cutaways/{root,transition,max_chord,span70,tip}_te.png
planform/{axes,structure}.png   (structure: web loci + chordwise rib lines)
full_blade/{elevation,closeups,closeups_te,cutaway,cutaway_te}.png
```

Regen also writes debug meshes at `figures/<key>/` (gitignored; same loft as the PNGs):

```
mids.obj / mids.vtk     full mapped loft (skins, webs, rib n-gons)
cutaway.obj / cutaway.vtk   same loft with suction skin omitted
samples.npz             station polylines, loft quads, all polygons + roles
```

![2D-F-LC-0R sections](figures/2D-F-LC-0R/cross_sections/sections.png)

![2D-F-LE-4R planform structure](figures/2D-F-LE-4R/planform/structure.png)

![2D-F-LE-10R close-ups](figures/2D-F-LE-10R/full_blade/closeups.png)

![2D-F-LE-10R close-ups TE](figures/2D-F-LE-10R/full_blade/closeups_te.png)

![2D-F-LE-10R elevation](figures/2D-F-LE-10R/full_blade/elevation.png)

![2D-F-LE-10R full blade without suction skin](figures/2D-F-LE-10R/full_blade/cutaway.png)

![3D-F-LC-10R full blade without suction skin](figures/3D-F-LC-10R/full_blade/cutaway.png)

![3D-F-LE-10R planform structure](figures/3D-F-LE-10R/planform/structure.png)

![4D-F-PA-10R planform structure](figures/4D-F-PA-10R/planform/structure.png)
