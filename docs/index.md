---
myst:
  html_meta:
    description: Signed-distance geometry for blade and wing sections, span lofts and conforming shell meshes.
---

# implicit-geometry-kernel

<div class="igk-hero">
<p class="lead">Geometry of blade and wing cross-sections and their span lofts:
aerofoils, planform, a catalogue of skin, spar-cap and shear-web layouts,
spanwise members, and a conforming shell mesh. Signed distance fields and
constructive solid geometry, with no materials, loads or solvers.</p>
</div>

![IEA 15 MW blade, true-scale close-ups from root to tip](../examples/blade/08_reference_blades/figures/iea15/full_blade/closeups.png)

```bash
pip install "implicit-geometry-kernel @ git+https://github.com/samhughes-fea/implicit-geometry-kernel@v0.9.5"
```

The documentation is four manuals, after the structure used by large open-source
engineering codes: **U** says how to use the kernel, **R** derives what it computes,
**V** shows the evidence that it computes it correctly, and **D** records how it is
developed. Pages are cited by key, for example [R1], and equations by number, (R1.4).

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} U · Use
:link: u/index
:link-type: doc

Installation, tutorials, how-to guides, the API and the layout catalogue.
:::

:::{grid-item-card} R · Reference theory
:link: r/index
:link-type: doc

Frames, aerofoils, signed distance fields, layouts, integration, loft and mesh,
derived equation by equation.
:::

:::{grid-item-card} V · Validation
:link: v/index
:link-type: doc

Analytic checks, convergence studies, golden numbers and reference blades.
:::

:::{grid-item-card} D · Development
:link: d/index
:link-type: doc

Architecture, workflow, testing, releases, decisions, roadmap and the bug log.
:::
::::

**Start here:** [getting started](tutorial.md) · [notation](r/R0_notation.md) ·
[changelog](project/changelog.md) ·
[GitHub](https://github.com/samhughes-fea/implicit-geometry-kernel)

```{toctree}
:caption: Manuals
:hidden:
:maxdepth: 2

u/index
r/index
v/index
d/index
```
