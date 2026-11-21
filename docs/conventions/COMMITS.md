# Commit messages

Every commit on `main` follows [Conventional Commits 1.0](https://www.conventionalcommits.org/en/v1.0.0/).
Pull requests are squash-merged, so the **PR title** becomes the commit on `main`;
CI checks the title. Commits on a branch should follow the same form.

## Form

```
<type>(<scope>): <summary>

<body: what changed and why, wrapped at 72 columns>

<footers>
```

- **Summary:** imperative, lower case after the colon, no full stop, under 72
  characters. Example: `fix(blade): classify a twisted flatback trailing edge as blunt`.
- **Body:** optional for small changes; say why, not how.
- **Footers:** `Closes #N` / `Refs #N`, `BREAKING CHANGE: ...`, co-author lines.

## Types

| Type | Use for | Changelog section |
|------|---------|-------------------|
| `feat` | new geometry capability, option or public name | Added |
| `fix` | wrong geometry, a crash, a failed mesh check | Fixed |
| `perf` | faster with the same results | Changed |
| `refactor` | restructuring with no behaviour change | Changed |
| `docs` | manual, docstrings, docs site, examples text | Documentation |
| `test` | tests only | — |
| `build` | packaging, `pyproject.toml`, dependencies | Changed |
| `ci` | workflows | — |
| `chore` | release prep, repository housekeeping | — |
| `revert` | reverting an earlier commit | as reverted |

## Scopes

Optional; use the label areas: `sdf`, `section`, `blade`, `structure`, `mesh`,
`docs-site`, `examples`, plus `release` for release-prep PRs and `deps` for
dependency bumps.

## Breaking changes

A removal or incompatible change of a public name adds `!` after the type or scope
(`feat(mesh)!: remove the flat mesh paths`) and a `BREAKING CHANGE:` footer naming
the replacement. Under [API stability](api/STABILITY.md) that only happens after a
release of deprecation warnings; the deprecation itself is a `feat` or `refactor`
with the `deprecation` label.

## Release PRs

`chore(release): prepare 0.9.N`. See [Tags and releases](RELEASING.md).
