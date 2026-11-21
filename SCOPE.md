# Scope

## What the kernel is

The source of truth for the **geometry** of blade and wing sections and their span
lofts. Downstream tools read this geometry and must not keep a parallel aerofoil,
planform or layout model.

## What it owns

- Aerofoils: coordinates, NACA and CST families, flatbacks, root blends.
- Planform: stations, chord, twist, pitch axis, offsets, span laws.
- The layout catalogue (`WX-Y-Z-NR`) and parametric layouts: skins, spar caps,
  shear webs, cores and inserts, with gauges in metres.
- Spanwise members: webs, ribs, stringers, bulkheads, the tip closure.
- Section cuts, the shell graph (mid-surface runs and their partitions), section
  geometric properties.
- The conforming analysis shell mesh: quads on skin and webs, triangles on ribs
  and the tip closure, with partition and thickness data on each element.
- Signed distance fields and CSG (`sdf/`), coordinate frames, integration, and
  readers that turn published formats (windIO, CSV, Selig) into geometry.

## What it never owns

- Materials, plies, laminate stiffness (ABD), section stiffness.
- Loads, solvers, element formulations, PCOMP cards or solver input decks.
- A Gmsh remesher or any other mesh driver beyond its own shell mesh.
- Data for specific confidential blades.

## Where it sits

Imports point down only:

```
implicit-geometry-kernel     geometry (this repository)
  -> laminate-constitutive-model   materials and laminates
    -> shell-fem                   3D shell finite elements
      -> anisotropic-section-model section stiffness
        -> 1d-fem                  beam solver
```

The kernel imports none of them.

## Deciding whether something belongs here

Ask in order:

1. Is it a property of shape, position or topology, measurable in metres or
   radians without knowing a material? If no, it goes downstream.
2. Would two different downstream tools both need it, and compute it the same
   way? If yes, it belongs here, so they share one answer.
3. Is a downstream adapter computing, converting or refusing geometry on top of
   the kernel today? Then it is a kernel gap; add it to the adapter gaps table in
   [the plan doc](docs/internal/REFACTOR_PLAN.md).
4. Does it make `src/` aware of one particular consumer? Then find the generic
   form or leave it downstream.

Data a consumer applies to geometry (for example a cap eccentricity used in a
constitutive law) is produced here as geometry and applied downstream.
