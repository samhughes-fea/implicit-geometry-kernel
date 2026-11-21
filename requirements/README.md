# Pinned CI dependencies

`ci.txt` pins every package CI installs (`[dev]` and `[docs]` extras), with
markers where Python 3.11 needs an older release. CI installs with
`pip install -c requirements/ci.txt -e ".[dev]"`, so a new upstream release cannot
break `main` overnight. The `core-install` job stays unpinned to catch that early.

Regenerate after changing dependencies in `pyproject.toml`:

```bash
uv pip compile pyproject.toml --extra dev --extra docs --universal --python-version 3.11 --no-header -o requirements/ci.txt
```

Dependabot proposes monthly bumps to this file.
