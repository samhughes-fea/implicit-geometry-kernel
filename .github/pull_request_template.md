<!-- Title: a Conventional Commit, e.g. "fix(blade): ..." (docs/conventions/COMMITS.md). It becomes the squash commit. -->

## What and why

Closes #

## Labels

One type, one or more areas, and any flag that applies (`deprecation`, `golden-change`, `consumer-gap`, `release`); see [LABELS](../docs/conventions/LABELS.md).

## Checks

- [ ] `pytest -q` passes; golden numbers unchanged, or `python scripts/snapshots.py --update` with the reason in `CHANGELOG.md`
- [ ] `ruff check src tests examples scripts`, `mypy src/implicit_geometry_kernel` and `python scripts/check_links.py` are clean
- [ ] `CHANGELOG.md` updated under the next `0.9.x` heading, line ending with this PR number
- [ ] `docs/internal/REFACTOR_PLAN.md` status updated
- [ ] No removed or changed public name without a release of `DeprecationWarning` first
- [ ] No private blade inputs (aerofoil coordinates, layups) and no consumer names in `src/` or `tests/`
- [ ] Regenerated example PNGs committed only where the example's subject changed
