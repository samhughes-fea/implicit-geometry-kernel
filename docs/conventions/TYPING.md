# Typing roadmap

mypy runs on all of `src/implicit_geometry_kernel` in CI. Strictness rises one
package at a time; a package moves up a level only when it passes at that level,
in the same PR that adds the annotations.

| Level | mypy settings | Packages |
|-------|---------------|----------|
| 3, strict | `disallow_untyped_defs`, `disallow_incomplete_defs`, `check_untyped_defs` | `structure`, `contract`, `frames`, `catalogue`, `extract` |
| 2, checked bodies | `check_untyped_defs`, `warn_unused_ignores` | `mesh` |
| 1, default | annotated functions only | `sdf`, `section`, `blade`, `integration`, top level |

`plotting` is excluded (Matplotlib types add noise and no safety).

## Next

1. `integration` to level 3 (nine signatures to annotate).
2. `sdf` to level 2, then 3: it is self-contained, so its types do not depend on
   the rest of the kernel.
3. `blade` and `section` to level 2. The loft mixins (`blade/loft/loft_*.py`)
   keep their override until the mixins become explicit protocols.
4. `mesh` to level 3.

## Rules

- New modules are written at the level of their package or higher.
- Use `numpy.typing.ArrayLike` for inputs that `np.asarray` accepts and
  `np.ndarray` for returns. `Any` is fine for duck-typed geometry arguments
  (`geo`, `cut`) until protocols exist.
- No blanket `# type: ignore`; give the code (`# type: ignore[attr-defined]`) and
  a reason.
