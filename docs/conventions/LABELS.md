# Labels

Labels are defined in [`.github/labels.yml`](../../.github/labels.yml) and synced to
GitHub from `main` by `label-sync.yml`. To add or change a label, edit that file in a
pull request; labels created by hand on GitHub are deleted on the next sync.

## Rules

- **One type label** on every issue and pull request: `bug`, `enhancement`,
  `documentation`, `maintenance` or `dependencies`. Issue templates and
  Dependabot set it; check it when a PR changes kind.
- **One or more area labels** (`area: sdf`, `area: section`, `area: blade`,
  `area: structure`, `area: mesh`, `area: docs-site`, `area: examples`), matching
  the package the change touches. Bug-log entries use the same areas.
- **Flags**, set whenever they apply:
  - `deprecation`: a public name is retired; the PR adds the warning and the
    [API stability](api/STABILITY.md) table row.
  - `golden-change`: golden numbers move on purpose; the changelog says why.
  - `consumer-gap`: closes or adds a row in the consumer adapter gaps table.
  - `release`: the release-prep PR (version bump and changelog heading; see
    [Tags and releases](RELEASING.md)).
  - `needs-decision`: waiting on the maintainer. Remove it once decided and
    record the decision in [docs/decisions/](../decisions/README.md) when it is
    a lasting one.
- **Closing labels**: `duplicate`, `wontfix` (link [SCOPE.md](../../SCOPE.md) when
  out of scope) or `question`.
- Agents apply labels when they open a PR and never remove `needs-decision`
  themselves.
