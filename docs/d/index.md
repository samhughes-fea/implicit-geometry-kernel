# D · Development

The D manual is for anyone who changes the kernel or wants to know how it is built,
tested, released and governed. It also holds the project's running record: the bug
and fix log, the feature development log and the tracker.

## Documents

| Key | Title | Status |
|---|---|---|
| D1 | Architecture and internals walkthrough | [Architecture](../architecture.md); traced call from `Station` to shell mesh planned |
| D2 | Contributing, workflow, commits and labels | [Contributing](../project/contributing.md), [shared standard](../conventions/standard/README.md), [conventions](../conventions/README.md), [rules](../project/rules.md) |
| D3 | Testing, CI jobs and scripts | Planned (see [tests/README](../../tests/README.md) and [.github/workflows](../../.github/workflows/)) |
| D4 | Releases and governance | [Releasing](../conventions/RELEASING.md), [security](../project/security.md), [code of conduct](../project/conduct.md) |
| D5 | Design decisions | [Architecture decision records](../decisions/README.md) |
| D6 | Roadmap | [Refactor plan](../internal/REFACTOR_PLAN.md) |
| D7 | Bug and fix log | [Section (SG)](../section/bug_log.md), [blade (BG)](../blade/bug_log.md), [frames (CF)](../frames/bug_log.md), [sdf (SD)](../sdf/bug_log.md); one index planned |
| D8 | Internal reports | [Internal notes](../internal/README.md) |
| D9 | Feature development log | Planned: one FD-nnn entry per feature |
| D10 | Tracker | Planned: generated from GitHub Issues at each build; today see [issues](https://github.com/samhughes-fea/implicit-geometry-kernel/issues) |

## Tracking

GitHub Issues is the one live tracker: it holds labels, assignees and `Closes #N`
links. The D manual shows its state (D10) and keeps the curated technical record:
D7 for defects (symptom, root cause, fix, test) and D9 for features (motivation,
design, R document changed, PRs, version, V check). Changelog lines link to the D7
and D9 entries they summarise.

```{toctree}
:maxdepth: 1
:caption: D documents

/architecture
/project/contributing
/conventions/README
/conventions/standard/README
/project/rules
/project/security
/project/conduct
/decisions/README
/internal/REFACTOR_PLAN
/section/bug_log
/blade/bug_log
/frames/bug_log
/sdf/bug_log
/internal/README
```
