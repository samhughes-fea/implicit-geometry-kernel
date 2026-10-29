# Feasibility and dropped webs

When resolve drops webs, steps the catalogue, or raises. Overview:
[`../README.md`](../../README.md). API table:
[`../api/README.md`](../../api/README.md#errors-and-feasibility).

↑ [`section_geometry`](../../README.md) · [`docs index`](../README.md)

## Table of contents

1. [Dropped webs vs infeasible](#dropped-webs-vs-infeasible)
2. [Sizing: freeze survivors](#sizing-freeze-survivors)
3. [Spanwise anchors](#spanwise-anchors)
4. [Zero-web keys `0A` / `0B`](#zero-web-keys)
5. [See also](#see-also)

---

<a id="dropped-webs-vs-infeasible"></a>
## Dropped webs vs infeasible

| Situation | Behaviour |
|-----------|-----------|
| Some anchored webs do not fit the local cavity | Those webs are **dropped**; listed on `StationLayoutResolution.dropped` with a `reason` (e.g. `no_inner_skin_intercept`) |
| Surviving web count maps to a lower catalogue rung | `layout.key` may differ from `nominal_layout.key`; `stepped_down` is `True` |
| Layout needs webs but **none** can be anchored | Raises `InfeasibleSectionGeometry` (not the same as “trial thicknesses too thick”) |
| Thicknesses exceed package clamps | Clamped before assembly; clamping alone does not raise |

Thickness trials that remain geometrically buildable after resolve are ordinary
failures only if a downstream consumer rejects them — geometry itself usually
still returns a `SectionGeometryResult`.

---

<a id="sizing-freeze-survivors"></a>
## Sizing: freeze survivors

Resolve once at minimum skin (caps/webs often zero), then pass
`resolution=` so later thickness trials keep the same web set:

Runnable: [`../examples/03_resolve_then_size.py`](../../examples/03_resolve_then_size.py).
Recipe: [`caller_patterns.md`](caller_patterns.md#2-resolve-once-then-rebuild-sizing-loop).

---

## Spanwise anchors

`web_x_positions_m` holds absolute offsets from the LE (typically fixed at the
root). On a short tip chord, outboard webs drop; they are **not** slid LE-ward.

Omit `web_x_positions_m` for the legacy `frac × local_chord` form — webs scale
with the station and are not dropped for taper.

Runnable comparison:
[`../examples/04_spanwise_web_anchors.py`](../../examples/04_spanwise_web_anchors.py)
(`5D-F` tip vs root-anchored positions).

---

<a id="zero-web-keys"></a>
## Zero-web keys `0A` / `0B`

| Key | Geometry | Shell-mesh handoff |
|-----|----------|--------------------|
| `0A` | Outer skin only | Intended (skin strips; no web junctions) |
| `0B` | Skin + single upper/lower spar-cap band | Intended (skin + cap strips; no web junctions) |

Both are valid catalogue keys for solids / properties **and** shell mesh.
Web feet and cap–web interiors are simply absent when `n_webs == 0`. Historical
“unsupported” wording described junction code that assumed webs — an
implementation gap, not the product rule. See
[`ADR_ZERO_WEB_SHELL_MESH.md`](../reference/ADR_ZERO_WEB_SHELL_MESH.md) and
[`geometry_to_mesh.md`](geometry_to_mesh.md).

---

## See also

- Layout cascade: [`../layout/`](../../layout/README.md)
- Taxonomy / step-down: [`../reference/system_type_xyz_taxonomy.md`](../reference/system_type_xyz_taxonomy.md)
- Quickstart sizing / spanwise: [`quickstart.md`](quickstart.md)
