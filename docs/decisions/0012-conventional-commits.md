# 0012. Commits follow Conventional Commits

**Decision.** Every commit on `main` follows Conventional Commits 1.0. Pull requests
are squash-merged, the PR title becomes the commit, and CI checks the title. Types,
scopes and the breaking-change marker are in [Commit messages](../conventions/standard/commits.md).

**Why.** History that says what kind of change each commit is, and where, can be
read, searched and turned into release notes without opening every diff.

**Follows.** Scopes match the label areas. A `!` marks a breaking change and only
follows a release of deprecation warnings.
