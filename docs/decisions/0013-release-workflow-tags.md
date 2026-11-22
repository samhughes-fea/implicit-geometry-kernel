# 0013. The release workflow makes tags

**Decision.** Merging a release PR (`chore(release): prepare 0.9.N`) into `main` is
the release. The `release` workflow then creates the annotated `v0.9.N` tag on that
merge commit and makes the GitHub release, using `GITHUB_TOKEN`. A person or agent
no longer creates tags by hand; a hand-made tag is still checked and released, but
it is the fallback. This supersedes "only the maintainer creates tags" in
[0008](0008-versioning.md) and [0011](0011-tags.md); every other rule in them stands.

**Why.** The maintainer dropped the maintainer-only rule on 2026-11-22. Agents
cannot push tags, and the hand step went wrong once: `v0.9.4` was made on the commit
before its release PR merged, so it carries 0.9.3 code. A workflow that reads the
version from the merged commit cannot tag the wrong commit.

**Follows.** The workflow refuses any version that is not `0.9.N` or `0.9.NrcM`, so
it never makes `v1.0.0` or `v0.10.0`. It does nothing when the tag already exists.
A release-PR merge that was not tagged is released by running the workflow by hand
with the commit SHA; the same hand run finishes a release whose tag was made but
whose GitHub release or upload failed. Each push to `main` is checked commit by
commit, so a release PR is found whether it is the head of the push or not, and
there is no workflow-wide concurrency group that could cancel a waiting release
run (22 Nov review). Rules and steps: [Tags and releases](../conventions/RELEASING.md).
