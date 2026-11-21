# Workflow: from issue to release

Every change can be traced from the release that shipped it back to the reason it
was made.

```
issue  ->  branch  ->  pull request  ->  squash commit on main  ->  release tag
 #N        any name    "Closes #N"       conventional message       v0.9.N
```

## 1. Issue

- Open an issue before non-trivial work: a bug, a feature, or a consumer gap. Use
  the templates; label it per [Labels](LABELS.md).
- Put it in the milestone of the release that should ship it (`0.9.N`). The
  maintainer creates one milestone per planned release.
- Trivial changes (a typo, a broken link) can skip the issue.

## 2. Pull request

- Title: a [Conventional Commit](COMMITS.md); it becomes the commit on `main`.
- Body: what and why, and `Closes #N` (or `Refs #N`) for each issue it answers.
- Labels and milestone match the issue.
- CI green, CODEOWNERS review, then squash-merge.

## 3. Records that point back

- **Changelog.** Each entry ends with the PR number, for example `(#81)`.
- **Bug logs.** A `docs/section/bug_log.md` or `docs/blade/bug_log.md` entry
  names its issue and fixing PR in the **Commit** line.
- **Plan doc.** A row in `docs/internal/REFACTOR_PLAN.md` that tracks the work
  links the issue (`#N`) so the table and GitHub agree.
- **Decisions.** A lasting choice settled on an issue gets a record in
  [docs/decisions/](../decisions/README.md) that links the issue.

## 4. Release

The release PR closes the milestone; the tag and GitHub release follow
[Tags and releases](RELEASING.md). The GitHub release notes are the changelog
section, so every line links to its PR.
