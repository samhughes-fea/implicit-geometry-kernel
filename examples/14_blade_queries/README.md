# 14 · Thickness field, queries, OML mesh and validate

Ply-drop ramps add geometry-only thickness T(z, s). `query` batches member,
thickness, OML distance and a local frame. `write_oml` writes a watertight
STL/OBJ. `validate()` names the station or member that failed.

```bash
python examples/14_blade_queries/blade_queries.py
```
