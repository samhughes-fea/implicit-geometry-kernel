# API stability

The kernel is pre-1.0 and stays on `0.9.x` while the public API settles. Within
that line, these rules hold.

## What is public

- Names in the `__all__` of the top-level package and of each subpackage
  (`blade`, `section`, `structure`, `sdf`, `mesh`, `extract`, `frames`,
  `catalogue`, `contract`, `integration`, `plotting`). These are the names the
  [API reference](../../api/index.md) documents.
- Their call signatures, return types and documented behaviour.

Everything else is internal: modules not re-exported by a package, names starting
with `_`, and attributes reached through a package `__getattr__` for tests. They
can change in any release. A monkeypatch of an internal name must target the
module that uses it, which can move.

## Changing a public name

1. Add the new name or behaviour. Adding is allowed in any release.
2. Keep the old name working for at least one released `0.9.x` with a
   `DeprecationWarning` that names the replacement and the release it goes in.
3. Remove it in a later `0.9.x` and list it under "Removed" in `CHANGELOG.md` and
   in [MIGRATING](../../MIGRATING.md).

Bug fixes that change results (for example a corrected trailing-edge
classification) are not deprecations; they are listed in the changelog with the
bug-log entry.

## Checks

- `tests/test_snapshots.py` compares the call signatures of the top-level
  `__all__`, `contract` and `structure` with `tests/golden/api.json`; a removed
  or changed signature fails it.
- Golden numbers (`tests/golden/sections.json`, `blade.json`) catch geometry that
  moves without a changelog entry.
- The docs build runs with `-W`, so a public name whose docstring does not render
  fails CI.

## Currently deprecated

| Name | Use instead | Goes in |
|------|-------------|---------|
| `mesh.shell_mesh` (module), `mesh.vtk`, `mesh.polygon`, `mesh.polygon_voronoi`, `mesh.delaunay` | `mesh.shell`, `mesh.export`, `mesh.fill` | 0.9.5 |
| Rare topology helpers imported from the package root | `implicit_geometry_kernel.blade` / `blade.topology` | announced in a later 0.9.x |
