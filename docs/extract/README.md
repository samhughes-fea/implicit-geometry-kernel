# extract

In-plane geometry extractors. This is a **subpackage** of
`implicit_geometry_kernel`, not a kernel-barrel export: import
`implicit_geometry_kernel.extract`, not `from implicit_geometry_kernel import ...`.

It returns **polylines**, not meshes. A 3-D FEM mesh stays downstream.

## Public names

| Symbol | Role |
| --- | --- |
| `extract_zero_loops` | Closed (or open) zero-level polylines of a 2-D SDF |
| `feature_curve` | Tagged feature curve from the 2.5-D field |
| `classify_trailing_edge` | Trailing-edge classification helper |

Callers speak `SDFNode` / the feature field, not section `Expr`. Station geometry
callers never need this package.

```python
from implicit_geometry_kernel.extract import extract_zero_loops, feature_curve
```
