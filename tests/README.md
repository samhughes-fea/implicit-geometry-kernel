# tests

```bash
pytest -q
```

↑ [repository](../README.md) · [docs](../docs/README.md)

Package folders mirror `src/implicit_geometry_kernel/`. Root-level tests cover loft
oracles, import layers, snapshots, and naming rules.

| Folder | Covers |
|---|---|
| [blade](blade/README.md) | Loft paths, moved modules |
| [catalogue](catalogue/README.md) | Tokens, cascade |
| [section](section/README.md) | Station worker, layouts, products |
| [structure](structure/README.md) | Assembler, contributions |
| [sdf](sdf/README.md) | SDF / CSG identity |
| [frames](frames/README.md) | Frame maps |
| [integration](integration/README.md) | Quadtree / shoelace |
| [contract](contract/README.md) | DTO shape |
| `golden/` | Number snapshots (`python scripts/snapshots.py`) |

Cross-cutting: `test_import_boundary.py`, `test_snapshots.py`,
`test_consumer_names.py`, `test_shell_graph.py`.

Examples are not pytest; CI runs `python scripts/run_examples.py`.
