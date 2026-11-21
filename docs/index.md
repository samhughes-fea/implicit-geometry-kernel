# implicit-geometry-kernel

This HTML site is `sphinx-build` over the **same** Markdown already in the
repository, plus an [API reference](api/index.md) generated from the docstrings
of each package's public names. GitHub still renders the Markdown in the code
tab; there is no `docs/manual/` copy.

Start at the [documentation hub](README.md) or the [glossary](glossary.md).

```{toctree}
:caption: Manual
:maxdepth: 1

README.md
glossary.md
sdf/README.md
section/README.md
blade/README.md
frames/README.md
structure/README.md
catalogue/README.md
contract/README.md
extract/README.md
integration/README.md
api/index.md
```

```{toctree}
:caption: Internal
:maxdepth: 1

conventions/README.md
internal/REFACTOR_PLAN.md
```
