# 0002. `sdf/` stays self-contained

**Decision.** Nothing in `sdf/` imports `section/` or `blade/`.

**Why.** The signed distance and CSG toolbox is useful on its own and could become
a general implicit-geometry library later.

**Follows.** Blade-specific shapes adapt into `sdf/` from the blade side
(for example `Airfoil.to_sdf`), never the other way.
