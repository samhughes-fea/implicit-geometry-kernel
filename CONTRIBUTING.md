# Contributing

Thanks for helping. Read [AGENTS.md](AGENTS.md) (the rules) and
[SCOPE.md](SCOPE.md) (what belongs here) first. The kernel owns blade and wing **geometry** only: aerofoils,
planform, the layout catalogue, gauges in metres, spanwise members, span topology
and the shell mesh. Materials, stiffness, loads and solvers belong in packages
built on top of it.

## Set up

```bash
git clone https://github.com/samhughes-fea/implicit-geometry-kernel
cd implicit-geometry-kernel
pip install -e ".[dev]"
pre-commit install          # optional: ruff and the link check on each commit
```

GEOS is needed for polygon offset and overlay; see [Install](README.md#install).

## Checks

Run these before opening a pull request. CI runs the same ones.

```bash
ruff check src tests examples scripts
mypy src/implicit_geometry_kernel
python scripts/check_links.py
pytest -q                               # about 8 minutes
python scripts/snapshots.py             # golden numbers and public-API snapshot
python scripts/run_examples.py          # optional: every example script
sphinx-build -W -b html docs docs/_build/html   # with the [docs] extra
```

## Rules

- **Golden numbers.** A refactor leaves `tests/golden/*.json` unchanged. A change
  that moves geometry on purpose runs `python scripts/snapshots.py --update` and
  says why in `CHANGELOG.md`.
- **Public API.** See [API stability](docs/conventions/api/STABILITY.md). Adding a
  public name is fine; removing or changing one needs a release of
  `DeprecationWarning` first.
- **Versions and tags.** Releases stay on `0.9.x`. Follow
  [Tags and releases](docs/conventions/RELEASING.md): the release PR bumps
  `__version__` and the changelog heading; the maintainer makes the tag.
- **Plan doc.** Update the status in `docs/internal/REFACTOR_PLAN.md` in the same
  PR as the work.
- **Generic code.** `src/` and `tests/` do not name downstream projects or special
  case them.
- **No private data.** Never commit confidential blade inputs (aerofoil
  coordinates, layups, planforms). Public reference blades are read from a
  separate checkout set by `REFERENCE_BLADES_DIR`, never copied in.
- **Figures.** Example PNGs are committed and shown in the example READMEs. Font
  rendering differs between machines, so commit a regenerated PNG only when the
  example's subject changed. Mesh exports (STL, OBJ, VTK) are not committed.
- **Docs.** The manual is the Markdown under `docs/`; the API reference is built
  from docstrings (NumPy style). Keep both building with `-W`.

## Bugs

Log geometry defects in `docs/section/bug_log.md` or `docs/blade/bug_log.md`
with symptom, root cause, fix and the test that checks it.
