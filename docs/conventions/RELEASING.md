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
- **Made by the workflow.** Merging the release PR is the release:
  `release.yml` creates the tag and the GitHub release. Nobody pushes a tag by
  hand unless the workflow cannot (see [Fallback](#fallback-a-tag-by-hand)).
  Contributors and agents prepare the release PR; whoever merges it releases.

## Release steps

1. A release PR, titled `chore(release): prepare 0.9.N` (the title is the
   squash commit, and the workflow looks for it), bumps `__version__` to
   `0.9.N`, renames the changelog section to `## 0.9.N (YYYY-MM-DD)`, updates
   [MIGRATING](../MIGRATING.md) when anything is deprecated or removed, updates
   `CITATION.cff` `version`, and records the release in the plan doc status
   table.
2. CI is green and the PR is squash-merged into `main`.
3. On that push to `main`, `release.yml` looks at every commit the push added
   to `main` (not only the head) for one titled `chore(release): prepare 0.9.N`:
   the squash subject, or the PR title on the first body line of a
   `Merge pull request #N` commit. A push that carries two release PRs fails;
   run the workflow by hand on each, oldest first. It then checks, against the
   release commit:
   - the subject's version equals `__version__`;
   - the version is `0.9.N` or `0.9.NrcM`; anything else, including `1.0.0` or
     `0.10.0`, fails the run and makes no tag;
   - the tag `v0.9.N` does not exist yet; if it already sits on that commit,
     the run stops and changes nothing, and if it sits on another commit, the
     run fails (tags never move; release the next patch);
   - `CHANGELOG.md` has a `## 0.9.N (YYYY-MM-DD)` heading (a release candidate
     needs only `## 0.9.N`).
4. It builds the sdist and wheel, then creates the annotated tag
   (`implicit-geometry-kernel 0.9.N`) on the merge commit and pushes it with
   `GITHUB_TOKEN`.
5. It makes the GitHub release, with the changelog section as its notes and the
   sdist and wheel attached (a release candidate is a pre-release), and uploads
   to PyPI only when `PYPI_PUBLISH` is `true`.
6. Consumers move their pins from the old tag to the new one in their own PRs.

Any other push to `main` stops after the first check; nothing is tagged.

A tag pushed with `GITHUB_TOKEN` does not start a new workflow run, so the same
run does the build and release; it does not wait for a tag push.

Runs do not share one workflow-wide concurrency group. GitHub keeps only one
pending run per group, so with one group a later push to `main` could cancel a
release run that was still waiting. The tag, GitHub release and PyPI jobs each
have a group per tag instead, so only two runs for the same release wait on
each other.

### A release PR merged without a tag

If a release PR merged before this workflow existed, or its run failed before
the tag was made, run the workflow by hand: **Actions > release > Run
workflow**, with the merge commit SHA. The run refuses a commit that is not on
`main` or that is not a release PR merge.

### A release that failed after its tag was made

If the GitHub release or the PyPI upload failed after the tag was pushed, either
use **Re-run failed jobs** on that run, or run the workflow by hand with the
merge commit SHA. A hand run that finds the tag already on that commit leaves
the tag alone and builds, makes or updates the GitHub release, and uploads to
PyPI (files already on PyPI are skipped). A tag on a different commit fails the
run.

### Fallback: a tag by hand

Only when the workflow cannot run (for example, Actions is down). Tag the
release PR's merge commit, never a later commit:
`git tag -a v0.9.N -m "implicit-geometry-kernel 0.9.N" <merge-sha> && git push origin v0.9.N`.
The tag push runs `release.yml`, which checks the same rules and makes the
release. Check that the target commit's `__version__` is `0.9.N`: the v0.9.4
tag was made by hand on the commit before the release PR merged and so carries
0.9.3 code.

The rule for `main` between releases is: `__version__` is the next release
number, and the changelog's top section collects the changes for it.

## Pinning

- Consumers pin a tag: `implicit-geometry-kernel @ git+https://github.com/samhughes-fea/implicit-geometry-kernel@v0.9.N`.
- A git submodule points at the tagged commit, and the consumer's changelog says
  which tag.
- Never pin a branch.

## Repository settings (maintainer)

- **Workflow permissions.** `release.yml` asks for `contents: write` only in the
  jobs that push the tag and make the release; the repository's default token
  can stay read-only.
- **Tag ruleset.** A tag ruleset on `v*` (Settings > Rules > Rulesets > New tag
  ruleset) with **Restrict updates** and **Restrict deletions** makes the
  "permanent" rule hold on GitHub as well as on paper. If it also has
  **Restrict creations**, add **GitHub Actions** to its bypass list, or the
  workflow cannot create the tag.
