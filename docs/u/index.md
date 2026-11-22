# U · Use

The U manual is for anyone who calls the kernel: how to install it, how to build a
blade and cut a section, and where every public name, code and data product is
listed. It states what to do; the reasons and the mathematics are in the
[R manual](../r/index.md), and the evidence that the results are correct is in the
[V manual](../v/index.md).

## Documents

| Key | Title | Status |
|---|---|---|
| U0 | Introduction: purpose, scope, what the kernel does not do, how to cite | [Documentation hub](../README.md) and [scope](../project/scope.md); to be rewritten |
| U1 | Installation: optional extras, platforms, system libraries | Planned (install notes are in the repository README) |
| U2 | Tutorials | [Getting started](../tutorial.md) |
| U3 | How-to guides | [Section quick start](../section/guides/quickstart.md), [anatomy on a figure](../section/guides/anatomy_on_a_figure.md), [caller patterns](../section/guides/caller_patterns.md), [geometry to mesh](../section/guides/geometry_to_mesh.md), [feasibility and drops](../section/guides/feasibility_and_drops.md) |
| U4 | Reference: API, layout catalogue, data contract, numerical defaults, errors, glossary | [API](../api/index.md), [catalogue](../catalogue/README.md), [contract](../contract/README.md), [glossary](../glossary.md); defaults and error index planned |
| U5 | Stability, deprecation, migration and changelog | [API stability](../conventions/api/STABILITY.md), [migrating](../MIGRATING.md), [changelog](../project/changelog.md) |
| U6 | Known limitations and performance | Planned |

```{toctree}
:maxdepth: 1
:caption: U documents

/README
/project/scope
/tutorial
/section/guides/quickstart
/section/guides/anatomy_on_a_figure
/section/guides/caller_patterns
/section/guides/geometry_to_mesh
/section/guides/feasibility_and_drops
/api/index
/catalogue/README
/contract/README
/glossary
/MIGRATING
/project/changelog
```
