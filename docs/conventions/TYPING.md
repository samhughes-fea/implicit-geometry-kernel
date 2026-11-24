# Typing roadmap

mypy runs on all of `src/implicit_geometry_kernel` in CI. Strictness rises one
package at a time; a package moves up a level only when it passes at that level,
in the same PR that adds the annotations.

| Level | mypy settings | Packages |
|-------|---------------|----------|
| 3, strict | `disallow_untyped_defs`, `disallow_incomplete_defs`, `check_untyped_defs` | `structure`, `contract`, `frames`, `catalogue`, `extract`, `integration` |
| 2, checked bodies | `check_untyped_defs`, `warn_unused_ignores` | `mesh`, `sdf`, `section`, `blade` (loft mixins keep their override) |
| 1, default | annotated functions only | `geos`, `profiling`, top level |

`plotting` is excluded (Matplotlib types add noise and no safety).

## Next

1. `sdf` to level 3: it is self-contained, so its types do not depend on the rest
   of the kernel (about 150 signatures).
2. `mesh` to level 3.
3. `section` and `blade` to level 3. The loft mixins (`blade/loft/loft_*.py`)
   keep their override until the mixins become explicit protocols.
4. `geos` to level 2 (two `has-type` errors in `geos/geom.py`).

## Rules

The shared [typing standard](standard/typing.md).
