# layout — layout identifiers and web cascade

↑ [`section_geometry`](../README.md)

Resolves a layout id such as `"2D-F"` (web count · spar-cap family ·
orientation) to a frozen `SystemLayoutSpec`, and owns the taper-time web-drop
cascade. Solid assembly is [`../sections/`](../sections/README.md). Decode:
[`../docs/reference/system_type_xyz_taxonomy.md`](../docs/reference/system_type_xyz_taxonomy.md).

Absolute web x-offsets [m from LE] may come from spanwise anchoring. Webs
whose local cavity is too shallow are **dropped**, not slid.
`effective_layout_for_web_count` maps survivors onto a catalogue rung;
`step_down_cascade` returns the ladder.

## Table of contents

1. [Organisation](#organisation)
2. [API](#api)
3. [Usage in the package](#usage-in-the-package)
4. [See also](#see-also)

---

## Organisation

| File | Role |
|------|------|
| `system_layout.py` | `SystemLayoutSpec`, registry, resolve, cascade |
| `structural.py` | Cap-family and fixed-cap-anchor parsers |
| `cap_placement.py` | `CapPlacementRule` — declarative cap policy |
| `section_layout.py` | Parametric `SectionLayout` facade |
| `section_axes.py` | Pitching-axis and max-thickness helpers |
| `thickness_limits.py` | Minimum cavity-depth hints |
| `types.py` | `SDFCallable` protocol |

## API

| Name | File | Role |
|------|------|------|
| `SystemLayoutSpec` | `system_layout.py` | Frozen layout config |
| `resolve_system_type` | `system_layout.py` | String → spec |
| `SYSTEM_TYPE_KEYS` | `system_layout.py` | Registered ids |
| `effective_layout_for_web_count` | `system_layout.py` | Surviving webs → catalogue rung |
| `step_down_cascade` | `system_layout.py` | Ordered drop ladder |
| `StructuralFamily` | `structural.py` | Cap family `"A"`…`"D"` |
| `CapPlacementRule` | `cap_placement.py` | Cap kind / anchor / half-width |
| `SectionLayout` | `section_layout.py` | Parametric web + orientation + caps |

### `SectionLayout` / `CapPlacementRule`

`SectionLayout.from_catalogue_key` builds a parametric view of a catalogue
entry: web chord fractions, orientation, and a `CapPlacementRule` derived from
the structural family (`none` / `fixed_band` / `per_web_pads` /
`continuous_box`). Zero-web family-B stations fold the catalogue half-band
fraction into `CapPlacementRule.half_width` before `build_spar_caps`.
Catalogue keys are presets; `multicell_from_section_layout` /
`multicell_from_layout(SectionLayout, …)` assemble without requiring a key.
`SectionGeometryRunner` / `build_section_geometry` /
`resolve_station_layout` accept the same `SectionLayout`: keyless drop is
``N`` surviving webs + `CapPlacementRule` (continuous box demotes when
``N < 2``); keyed layouts still use the catalogue cascade.
`SystemLayoutSpec` remains the registry source of truth for named presets.

### `SystemLayoutSpec` fields

| Field | Meaning |
|-------|---------|
| `key` | Catalogue string (e.g. `"2D-F"`) |
| `n_webs`, `n_cells` | Web and cell counts |
| `has_spar_caps` | Whether caps are present |
| `web_orientation` | `"chord_normal"` or `"flapwise"` |
| `web_positioning` | How default chord fractions are chosen |
| `web_chord_fracs` | Default web stations as fractions of chord from LE |
| `spar_cap_positioning` | Cap anchor policy |
| `geometry_mode` | Catalogue metadata; `"zero_web_catalogue"` on `0A`/`0B` (host is still `MultiCellSection`) |
| `structural_family` | Cap family `"A"`…`"D"` |
| `description` | Human-readable note |

## Usage in the package

```mermaid
flowchart LR
  Key["2D-F"] --> Resolve[resolve_system_type]
  Resolve --> Spec[SystemLayoutSpec]
  Spec --> Station[resolve_station_layout]
  Station --> Sections[sections]
```

Station resolve is orchestrated from [`../api/`](../api/README.md)
(`resolve_station_layout`).

**Visual SoT for root-anchored vs station-relative webs:** the demo
[`../visualisation/demos/spanwise_yz_plan_flyover.py`](../visualisation/demos/spanwise_yz_plan_flyover.py)
animates absolute LE tracks along the planform and shows webs **dropped** when
the local cavity can no longer hold them (not slid). BSP L2 stamps
`web_x_positions_m`; this package owns the drop / catalogue step-down.
`section_schedule` does **not** own this policy (thickness presence is a
different spanwise product — see BSP
[`glossary`](../../../../blade_structure_precompute/module/section_selection/docs/glossary.md)).

## See also

- [`../docs/reference/system_type_xyz_taxonomy.md`](../docs/reference/system_type_xyz_taxonomy.md)
- [`../sections/`](../sections/README.md)
- [`../visualisation/demos/spanwise_yz_plan_flyover.py`](../visualisation/demos/spanwise_yz_plan_flyover.py)
- BSP [`package_boundaries.md`](../../../../blade_structure_precompute/module/section_selection/docs/package_boundaries.md) §3.2 / §3.10
