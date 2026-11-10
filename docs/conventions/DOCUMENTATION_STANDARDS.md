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

The Sphinx + MyST site (`docs/conf.py`, extra `[docs]`) consumes these Markdown
files. It must not become a second manual: no `sphinx.ext.autodoc`, no
`docs/manual/` copy. GitHub still renders the same files in the code tab.

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

Package logs stay at that level (not nested further). They are internal:

| Defect lives in | Log | ID |
|---|---|---|
| Span loft, mapped strips, station CSG viz | [`docs/blade/bug_log.md`](../blade/bug_log.md) | `BG-nnn` |
| Section geometry, midlines, junctions | [`docs/section/bug_log.md`](../section/bug_log.md) | `SG-nnn` |
| Frames / transforms | [`docs/frames/bug_log.md`](../frames/bug_log.md) | `CF-nnn` |

A confirmed defect is not closed until the matching log has an index row and a
body. Chat, `CHANGELOG.md`, and a guide are not a substitute.

Write the entry **when the bug is found**, not after the patch. Update the same
entry with evidence from the fix, then mark it resolved. Do not delete closed
entries; if a first diagnosis was wrong, strike or narrate it.

| When | Required fields |
|---|---|
| Found | Next ID; index row **Open**; **Symptom**; **Evidence** (figure, station `z`, test, metric) |
| Fix in hand | **Root cause**; **Fix** (behaviour and files); tests or figures named |
| Resolved | Status **Fixed**; **Verification** (the test or regenerated figure that proves it) |

A follow-on of the same mechanism updates that ID. A distinct symptom gets a new
ID (do not bury 5D CSG hold under an IEA wrap rewrite). Do not log an intended
design change or a figure regen that was not a defect.

Agents bound by this file must open or update the log in the same change set as
the code fix.

## Consumer names

Live `src/` and `tests/` must not name consuming packages
(`tests/test_consumer_names.py`). `docs/FORK.md` and scripts may.
