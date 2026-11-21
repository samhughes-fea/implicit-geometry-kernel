# Tags and releases

A tag is a promise that consumers can pin. These rules keep every tag meaningful
and permanent. `release.yml` enforces the checkable ones.

## What a tag is

- **Format:** `v0.9.N` for a release, `v0.9.NrcM` for a release candidate
  ([0008](../decisions/0008-versioning.md): no 0.10.0, no 1.0.0 until the
  maintainer decides). Nothing else starting with `v`.
- **Only releases are tagged.** No working tags (such as `K9`) for in-progress
  work. To pin unreleased work, a consumer pins a commit SHA on `main`, and moves
  to the next release tag as soon as it exists.
- **Annotated**, with the message `implicit-geometry-kernel 0.9.N`.
- **On `main`.** The tagged commit is the merge of the release PR, never a
  branch commit.
- **Matches the code.** The tag equals `v` + `__version__`, and `CHANGELOG.md` has
  a `## 0.9.N (date)` heading.
- **Permanent.** A pushed tag is never moved, re-pointed or deleted. A mistake is
  fixed by the next patch release (`0.9.N+1`), and the bad release is marked in
  the changelog.
- **Maintainer only.** Only the maintainer creates tags and GitHub releases.
  Contributors and agents prepare the release PR and stop there.

## Release steps

1. A release PR (`chore(release): prepare 0.9.N`) bumps `__version__` to `0.9.N`, renames the changelog section
   to `## 0.9.N (date)`, updates [MIGRATING](../MIGRATING.md) when anything is
   deprecated or removed, updates `CITATION.cff` `version`, and records the release
   in the plan doc status table.
2. CI is green and the PR is merged into `main`.
3. The maintainer tags the merge commit:
   `git tag -a v0.9.N -m "implicit-geometry-kernel 0.9.N" && git push origin v0.9.N`,
   or creates the release on GitHub's release page with a new `v0.9.N` tag on
   `main`.
4. `release.yml` checks the tag name, that the commit is on `main`, the
   changelog heading and `__version__`, then builds the sdist and wheel
   (and uploads to PyPI only when `PYPI_PUBLISH` is `true`).
5. The GitHub release notes are the changelog section.
6. Consumers move their pins from the old tag to the new one in their own PRs.

The rule for `main` between releases is: `__version__` is the next release
number, and the changelog's top section collects the changes for it.

## Pinning

- Consumers pin a tag: `implicit-geometry-kernel @ git+https://github.com/samhughes-fea/implicit-geometry-kernel@v0.9.N`.
- A git submodule points at the tagged commit, and the consumer's changelog says
  which tag.
- Never pin a branch.

## Repository settings (maintainer)

A tag ruleset on `v*` (Settings > Rules > Rulesets > New tag ruleset) with
**Restrict creations**, **Restrict updates** and **Restrict deletions**, and the
maintainer as the only bypass, makes the "maintainer only" and "permanent" rules
hold on GitHub as well as on paper.
