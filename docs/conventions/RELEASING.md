# Kernel releases

**Audience: internal.** The release rules, steps and pinning are the shared
[releasing standard](standard/releasing.md); `release.yml` is synced from it and
reads this repository's values from `.github/repo.yml`. This page keeps only what
is particular to the kernel.

## Version line

Every release is `0.9.N` (`version-line: "0.9"`,
[decision 0008](../decisions/0008-versioning.md)). The workflow refuses
`0.10.0` and `1.0.0`. Release PRs are titled `chore(release): prepare 0.9.N`.

## Changelog headings

From 0.9.6 the changelog uses the shared form: the top section is
`## [Unreleased]`, and the release PR renames it `## [0.9.N] - YYYY-MM-DD`.
Headings up to 0.9.5 keep their old `## 0.9.N (YYYY-MM-DD)` form.

## History

- `v0.9.4` was tagged by hand on the commit before its release PR merged, so it
  carries 0.9.3 code. Tags never move, so it stays; pin `v0.9.5` or later.
- `v0.9.5` was the first release tagged by the workflow
  ([decision 0013](../decisions/0013-release-workflow-tags.md)).

## Pinning the kernel

`implicit-geometry-kernel @ git+https://github.com/samhughes-fea/implicit-geometry-kernel@v0.9.N`.
A git submodule points at the tagged commit.
