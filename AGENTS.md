# Rules for agents and contributors

Read this before changing the repository. It applies to people and to coding
agents alike. [SCOPE.md](SCOPE.md) says what belongs here;
[ARCHITECTURE.md](ARCHITECTURE.md) says where things live;
[CONTRIBUTING.md](CONTRIBUTING.md) lists the checks; past choices are in
[docs/decisions/](docs/decisions/README.md).

## Hard rules

1. **Geometry only.** No materials, plies, stiffness, loads, solvers or a Gmsh
   driver in `src/`. If a feature needs one of those, it belongs downstream.
2. **Generic source.** `src/` and `tests/` never name a downstream project or
   special-case one (`tests/test_consumer_names.py` checks this). Public docs may
   name public sibling packages; nothing anywhere names a private project.
3. **No private data.** Never commit confidential blade inputs: aerofoil
   coordinates, layups, planforms or material cards from non-public blades.
   Reference blades are read from a separate checkout (`REFERENCE_BLADES_DIR`, or a
   git-ignored `reference_blades/`); tests and examples skip when it is absent.
   Figures and summary results of public reference blades are fine.
4. **Versions stay on 0.9.x.** Never set `__version__` to 1.0.0 or plan a 0.10.0.
   Bump the patch version and the `CHANGELOG.md` heading in the PR that prepares a
   release.
5. **Tags follow [the release rules](docs/conventions/RELEASING.md).** The
   `release` workflow creates tags and GitHub releases when a release PR
   (`chore(release): prepare 0.9.N`) is merged; nobody pushes a tag by hand
   unless the workflow cannot. A tag is `v0.9.N` (or `v0.9.NrcM`) on a `main`
   merge commit, equals `__version__`, has a dated changelog heading, and is
   never moved or deleted. No working tags: pin unreleased work by commit SHA.
   Agents prepare the release PR; merging it is the release.
6. **Deprecate before removing.** A public name (see
   [API stability](docs/conventions/api/STABILITY.md)) keeps working with a
   `DeprecationWarning` for at least one released 0.9.x before it goes.
7. **Golden numbers move only on purpose.** `python scripts/snapshots.py` must
   pass; a deliberate geometry change updates the snapshots and says why in the
   changelog.
8. **Plan doc in every PR.** Update the status in
   [docs/internal/REFACTOR_PLAN.md](docs/internal/REFACTOR_PLAN.md) in the same PR
   as the work, so anyone can resume from the repository alone.
9. **Fix the kernel, not the consumer.** Geometry a consumer has to compute,
   convert or refuse on top of the kernel is a kernel to-do: add it to the
   consumer adapter gaps table in the plan doc.

## Working conventions

- Follow [the workflow](docs/conventions/WORKFLOW.md): an issue for non-trivial
  work, `Closes #N` in the PR, the PR number on each changelog line, issue and PR
  in bug-log entries.
- Commit messages and PR titles follow [Conventional Commits](docs/conventions/COMMITS.md)
  (`fix(blade): ...`, `feat(mesh)!: ...`); PRs are squash-merged on their title.
- Label every issue and PR per [the label rules](docs/conventions/LABELS.md):
  one type, its areas, and any review flag. Never remove `needs-decision`.
- Work through pull requests; `main` is protected. Keep one concern per commit.
- Run the checks in [CONTRIBUTING.md](CONTRIBUTING.md) before pushing. After a
  failure, re-run only what failed while fixing it.
- `sdf/` stays self-contained: nothing in it imports `section/` or `blade/`.
- Example PNGs are committed and embedded in each example README. Commit a
  regenerated PNG only when the example's subject changed; font rendering differs
  between machines. Never commit mesh exports.
- Log a geometry defect in `docs/section/bug_log.md` or `docs/blade/bug_log.md`
  (symptom, root cause, fix, test) in the PR that fixes it.
- Docstrings are NumPy style and must build with `sphinx-build -W`.
- Write plainly: short sentences, units stated, no marketing words.
