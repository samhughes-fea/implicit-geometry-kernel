# Conventions

**Audience: internal.** How work is done in this repository.

↑ [docs](../README.md)

The shared development standard, synced into [standard/](standard/README.md),
holds the rules every tool repository follows: workflow, commits, pull requests,
labels, releasing, versioning, CI, testing, documentation, API stability,
typing, architecture, data and security, GitHub settings and planning. Do not
edit those files here.

The pages below are the kernel's own, and add to the standard.

| Doc | Role |
|---|---|
| [RELEASING.md](RELEASING.md) | Version line 0.9.x, changelog headings, release history |
| [api/STABILITY.md](api/STABILITY.md) | What is public in the kernel; currently deprecated names |
| [api/API_STANDARDS.md](api/API_STANDARDS.md) | Worker runners (`SectionGeometry.execute`) |
| [TYPING.md](TYPING.md) | Which package is at which mypy level, and what is next |
| [MODULE_ARCHITECTURE_STANDARDS.md](MODULE_ARCHITECTURE_STANDARDS.md) | The kernel's import layers, extras, dependencies |
| [DOCUMENTATION_STANDARDS.md](DOCUMENTATION_STANDARDS.md) | The kernel's doc layers, page types, bug-log IDs |

Commit scopes and area labels: [.github/repo.yml](../../.github/repo.yml) and
[.github/labels.yml](../../.github/labels.yml).
