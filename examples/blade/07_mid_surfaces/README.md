# 07 · Mid-surface patches

Lofted wall mid-surface patches as OBJ next to the OML STL.

```bash
python examples/blade/07_mid_surfaces/mid_surfaces.py
```

`BladeGeometry.write_oml` writes the outer mould line. `write_mid_surfaces` writes
one patch per loftable `shell_graph` run into `figures/` — the conformal **skin mid**,
not a stack-NA jump at cap ends. Optional extra `[cad]` exports STEP from `blade.cad`
only.

↑ [blade examples](../README.md)
