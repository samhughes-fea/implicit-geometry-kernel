# Extract

↑ [docs](../README.md) · [src](../../src/implicit_geometry_kernel/extract/README.md) ·
[tests](../../tests/README.md)

In-plane geometry extractors. Import `implicit_geometry_kernel.extract`; the
symbols are not on the kernel barrel. The package returns **polylines**, not
meshes. A three-dimensional finite-element mesh remains downstream.

## Zero set

For a planar field \(d(\mathbf{x})\), a zero loop is a polyline approximating

$$
\{\mathbf{x}\in\mathbb{R}^2: d(\mathbf{x})=0\}.
\tag{1}
$$

Function: `extract_zero_loops`. Tagged curves from the two-and-a-half-dimensional
feature field use `feature_curve`. `classify_trailing_edge` is a trailing-edge
helper. Callers speak `SDFNode` and the feature field, not section `Expr`.
In-plane Newton projection is `GeometrySolver` in `sdf.solve`; this package
traces polylines, it is not a second solver. Station geometry callers do not
need this package. Not on the kernel barrel; see [glossary](../glossary.md).

```python
from implicit_geometry_kernel.extract import extract_zero_loops, feature_curve
```
