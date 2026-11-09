# implicit_geometry_kernel

Installable package root. Public loft API is re-exported here; one-station geometry
is `implicit_geometry_kernel.section`.

↑ [repository](../../../README.md) · [docs](../../../docs/README.md) ·
[glossary](../../../docs/glossary.md)

| Subpackage | Docs | Examples | Tests |
|---|---|---|---|
| [blade](blade/README.md) | [docs/blade](../../../docs/blade/README.md) | [examples/blade](../../../examples/blade/) | [tests/blade](../../../tests/blade/) |
| [catalogue](catalogue/README.md) | [docs/catalogue](../../../docs/catalogue/README.md) | via section / blade concepts | [tests/catalogue](../../../tests/catalogue/) |
| [section](section/README.md) | [docs/section](../../../docs/section/README.md) | [examples/section](../../../examples/section/) | [tests/section](../../../tests/section/) |
| [structure](structure/README.md) | [docs/structure](../../../docs/structure/README.md) | [examples/structure](../../../examples/structure/) | [tests/structure](../../../tests/structure/) |
| [sdf](sdf/README.md) | [docs/sdf](../../../docs/sdf/README.md) | [examples/sdf](../../../examples/sdf/) | [tests/sdf](../../../tests/sdf/) |
| [frames](frames/README.md) | [docs/frames](../../../docs/frames/README.md) | [examples/frames](../../../examples/frames/) | [tests/frames](../../../tests/frames/) |
| [integration](integration/README.md) | [docs/integration](../../../docs/integration/README.md) | via section / sdf | [tests/integration](../../../tests/integration/) |
| [contract](contract/README.md) | [docs/contract](../../../docs/contract/README.md) | teaching handoff | [tests/contract](../../../tests/contract/) |
| [extract](extract/README.md) | [docs/extract](../../../docs/extract/README.md) | — | [tests/extract](../../../tests/extract/) |

Optional extras (not geometry layers; geometry must not import them):

| Extra | Docs | Notes |
|---|---|---|
| [plotting](plotting/README.md) (`[plots]`) | [visualisation](../../../docs/section/visualisation/README.md) | matplotlib |
| `mesh/` | polygonal fill; Gmsh in `mesh/gmsh/` (`[mesh]`) | geometry must not import this package |
| [blade/cad.py](blade/cad.py) (`[cad]`) | STEP writer | lazy cadquery |
| `profiling.py` | timers | `SECTION_GEOMETRY_PROFILE=1`; not on root `__all__` |
