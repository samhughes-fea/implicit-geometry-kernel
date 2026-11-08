# Documentation standards

Layered `README.md` files are the on-ramp. Guides, reference, and reports sit under
`docs/<package>/` and are linked from that package README.

↑ [conventions](README.md)

**Audience: internal** (this file is the house style for writers). The pages it
governs fall into two audiences below.

## Audiences

Keep one documentation tree. Do not duplicate a `docs/manual/` copy.

| Audience | Who | Location | Voice |
|---|---|---|---|
| **Manual** | Library users and developers of consuming tools | Package READMEs, guides, theory/reference, repository README, example indexes | Academic English; schematics; numbered equations |
| **Internal** | Maintainers of this repository | `docs/conventions/`, `bug_log.md`, `reports/`, ADRs, `REFACTOR_PLAN.md`, `FORK.md`, deep `docs/section/{assembly,cells,junctions,…}` pipeline pages | Terse; one-line **Audience: internal** banner |

`docs/README.md` owns two maps: Manual (capability table) and Internal (conventions, bug logs, reports).

`src/**/README.md` remain ten to thirty lines: purpose, key imports, links. They are not essays.

## Page types

One tree, four kinds of page. Do not duplicate a type across files.

| Kind | What it is | Where |
|---|---|---|
| **Overview** | What the package is, a mermaid, numbered equations, public attributes | `docs/<pkg>/README.md` |
| **How-to** | A task with a runnable twin | `docs/section/guides/`, numbered `examples/` |
| **Reference** | Look-up: identifiers, field tables, frames, names | taxonomy, [glossary](../glossary.md), frame cheat, package Results |
| **Internal** | Maintainer notes | bug logs, ADRs, reports, pipeline READMEs (banner required) |

**Completeness:** every name in `src/implicit_geometry_kernel/__init__.py` `__all__` appears in [docs/glossary.md](../glossary.md) with a link to one home page. Do not paste every dataclass into every README.

A later Sphinx, autodoc, or MkDocs site must consume these Markdown files. It must not become a second manual.

## Manual writing

### English

Write complete sentences. Define each symbol at first use. Prefer British academic
register (aerofoil, centreline, labelled). Chatty asides are allowed only when they
state a definition immediately afterwards.

### Schematics

Use mermaid for frame maps, data flow, catalogue decode, and loft versus station.
One diagram per conceptual hop. Node IDs are camelCase with no spaces.

Package tree, data flow, and example map live in `docs/README.md` and
`examples/README.md`. Section pipeline diagrams live under `docs/section/`.

### Equations

Use display math with tags. Follow each display with a sentence that names every
symbol. Place the equation next to the function that implements it.

```tex
$$
\mathbf{x}_B = \mathbf{R}(\theta)\, \mathbf{x}_S
\tag{1}
$$
```

where \(\theta\) is geometric twist about the span axis, \(\mathbf{x}_S\) is a
point in the zero-twist station frame, and \(\mathbf{R}(\theta)\) is the
corresponding plane rotation.

Do not rewrite historical reports to this standard.

## Layers

Every major directory has a `README.md` that links **up** (parent), **down**
(children), and **sideways** (matching `docs/`, `src/`, `tests/`, `examples/`).

| Layer | Location |
|---|---|
| Product | repository `README.md` |
| Doc hub | `docs/README.md` |
| Package | `docs/<pkg>/README.md` |
| Module | `docs/section/<module>/README.md` (pipeline pages; internal unless listed as a guide) |
| Code | `src/implicit_geometry_kernel/<pkg>/README.md` |
| Tests | `tests/<pkg>/README.md` |
| Examples | `examples/<pkg>/README.md` and `examples/<pkg>/NN_*/README.md` |

## Capability map

`docs/README.md` owns the capability table (API → products → docs → examples → tests).
Other READMEs repeat only the rows they own.

## Runnable twins

Guides that say “run this” must name a script that exists. Do not cite deleted
handoff or consumer-mesh drivers. The kernel stops at `contract` types; a
three-dimensional mesh is downstream.

## Bug logs

`docs/section/bug_log.md` and `docs/blade/bug_log.md` stay at that level (not nested
further). They are internal.

## Consumer names

Live `src/` and `tests/` must not name consuming packages
(`tests/test_consumer_names.py`). `docs/FORK.md` and scripts may.
