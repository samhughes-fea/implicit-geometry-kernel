# Rules for agents and contributors

<!-- guidelines:start (shared development standard 0.1.1; edit it there, not here) -->
## Shared rules

These rules hold in every tool repository; the full standard is in
[docs/conventions/standard/](docs/conventions/standard/README.md). The rules
under "Hard rules" below this block are this repository's own.

1. **Read before changing.** This file, `SCOPE.md`, then the plan doc named in
   `.github/repo.yml`.
2. **Work through pull requests.** `main` is protected. One concern per pull
   request; the title is a Conventional Commit
   ([commits](docs/conventions/standard/commits.md)) and becomes the squash commit.
3. **Keep the plan doc current** in the same pull request as the work.
4. **No private data.** Never commit confidential inputs or name a private
   project ([data and security](docs/conventions/standard/data-and-security.md)).
5. **Generic source.** `src/` and `tests/` never name or special-case a
   downstream project.
6. **Deprecate before removing** ([API stability](docs/conventions/standard/api-stability.md)).
7. **Golden numbers move only on purpose**, with the `golden-change` label and
   the reason in the changelog.
8. **Never remove `needs-decision`**; record a lasting decision in `docs/decisions/`.
9. **Releases follow the workflow.** Prepare the release pull request; never push
   or move a tag ([releasing](docs/conventions/standard/releasing.md)).
10. **Label everything**: one type, its areas, any flag
    ([labels](docs/conventions/standard/labels.md)).
11. **Example figures are committed** and embedded in the example README; commit a
    regenerated PNG only when the example's subject changed.
12. **Re-run only what failed** while fixing a failing check.
<!-- guidelines:end -->

Read this before changing the repository. It applies to people and to coding
agents alike. [SCOPE.md](SCOPE.md) says what belongs here;
[ARCHITECTURE.md](ARCHITECTURE.md) says where things live;
[CONTRIBUTING.md](CONTRIBUTING.md) lists the checks; past choices are in
[docs/decisions/](docs/decisions/README.md).

## Hard rules

These are this repository's own rules, numbered on from the shared ones above.

13. **Geometry only.** No materials, plies, stiffness, loads, solvers or a Gmsh
    driver in `src/`. If a feature needs one of those, it belongs downstream
    ([SCOPE.md](SCOPE.md)).
14. **Consumer names.** `tests/test_consumer_names.py` checks that `src/` and
    `tests/` name no downstream project. Public docs may name public sibling
    packages; nothing anywhere names a private project.
15. **Reference blades** are read from a separate checkout (`REFERENCE_BLADES_DIR`,
    or a git-ignored `reference_blades/`); tests and examples skip when it is
    absent. Figures and summary results of public reference blades are fine.
16. **Versions stay on 0.9.x** (`version-line: "0.9"` in `.github/repo.yml`,
    [decision 0008](docs/decisions/0008-versioning.md)). Never set `__version__`
    to 1.0.0 or plan a 0.10.0. Release PRs are `chore(release): prepare 0.9.N`
    ([kernel release notes](docs/conventions/RELEASING.md)).
17. **Golden numbers.** `python scripts/snapshots.py` must pass; a deliberate
    geometry change runs it with `--update` and says why in the changelog.
18. **The plan doc** is [docs/internal/REFACTOR_PLAN.md](docs/internal/REFACTOR_PLAN.md).
19. **Fix the kernel, not the consumer.** Geometry a consumer has to compute,
    convert or refuse on top of the kernel is a kernel to-do: add it to the
    consumer adapter gaps table in the plan doc, with the `consumer-gap` label.

## Working conventions

- Commit scopes and `area:` labels are in [.github/repo.yml](.github/repo.yml) and
  [.github/labels.yml](.github/labels.yml).
- Run the checks in [CONTRIBUTING.md](CONTRIBUTING.md) before pushing.
- `sdf/` stays self-contained: nothing in it imports `section/` or `blade/`.
- Never commit mesh exports (STL, OBJ, VTK).
- Log a geometry defect in `docs/section/bug_log.md`, `docs/blade/bug_log.md` or
  `docs/frames/bug_log.md` (symptom, root cause, fix, test) in the PR that fixes it.
- Docstrings are NumPy style and must build with `sphinx-build -W`.
