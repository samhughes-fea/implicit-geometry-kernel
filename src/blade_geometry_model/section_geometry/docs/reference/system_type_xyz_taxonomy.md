# Layout ids — how `"2D-F"` is decoded

Layout id (SystemType) pattern:

```text
{X}{Y}-{Z}     e.g.  2D-F
 │ │  └─ web orientation
 │ └─── spar-cap family Y
 └───── web count X
```

| Part | Meaning | Examples |
|------|---------|---------|
| **X** | Shear-web count | `0`, `1`, `2`, … |
| **Y** | Spar-cap family (`StructuralFamily`) | `A` none · `B` one band · `C` per-web bands · `D` continuous box |
| **Z** | Web orientation | `CN` chord-normal · `F` flapwise (omit when X = 0) |

Display form: `SystemType{X}{Y}-{Z}`. `resolve_system_type("2D-F")` → frozen
`SystemLayoutSpec`.

↑ [`section_geometry`](../../README.md) ·
[`../layout/README.md`](../../layout/README.md) ·
[`structural_systems_gallery.md`](structural_systems_gallery.md)

## Table of contents

1. [Y — spar-cap families](#y--spar-cap-families)
2. [Z — web orientation](#z--web-orientation)
3. [Layout registry (`resolve_system_type`)](#layout-registry-resolve_system_type)
4. [Spanwise step-down cascade](#spanwise-step-down-cascade)
5. [Reference axes and kwargs (advanced)](#reference-axes-and-kwargs-advanced)
6. [Downstream handoff](#downstream-handoff)
7. [Caveats (v1)](#caveats-v1)

---

## Y — spar-cap families

| Y | Description |
|---|-------------|
| **A** | No spar caps: skin, webs, cores (and optional inserts) only. |
| **B** | **Single** fixed upper/lower cap pair on **one** chordwise band. Anchor: **pitching** (default ~⅓ chord from LE) or **max_thickness**. |
| **C** | **Multiple** discrete upper/lower cap patches — one band per web — merged as `spar_cap_suction_surface` / `spar_cap_pressure_surface`. |
| **D** | **Continuous box** spar spanning leftmost to rightmost web (needs **X ≥ 2**; no `1D-*`). |

---

## Z — web orientation

- **CN** (chord-normal): chord-normal in the **untwisted** section.
- **F** (flapwise): vertical in the **twisted** blade frame.

Compare
[`../visualisation/output/structural_family_grids/2D-CN.png`](../../visualisation/output/structural_family_grids/2D-CN.png)
vs
[`../visualisation/output/structural_family_grids/2D-F.png`](../../visualisation/output/structural_family_grids/2D-F.png).

### Implementation note (S / B frames)

- **CN**: Web mid-surface runs chord-normal in the chord (**S**) frame before
  section twist.
- **F**: Chord-frame `ShearWeb` uses tilted S-frame endpoints so that, after
  the section’s global `+twist`, the axis is B-frame vertical and meets the
  same rotated inner-skin clip as skin and caps (`MultiCellSection` flapwise
  branch — not a `ShearWeb` counter-rotation).

Spar-cap clipping follows the same chord-vertical vs parallel-strip logic as
webs.

---

## Layout registry (`resolve_system_type`)

Source of truth:
[`../layout/system_layout.py`](../../layout/system_layout.py)
(`SYSTEM_TYPE_KEYS`, `_LAYOUT_REGISTRY`). The middle letter of each key is
taxonomy **Y** and matches `structural_family` on the spec.

- **0-web:** only **`0A`** and **`0B`** (**X = 0** ⇒ no **Z** suffix).
  Both build `MultiCellSection` with empty webs (family A / B). Catalogue
  field `geometry_mode="zero_web_catalogue"` is metadata only — it does not
  select a separate host class.
- **No `0C` / `0D`:** **Y = C** and **Y = D** need webs (`X ≥ 1`).
- **Multicell:** for **X = 1..5**, **Y = A, B, C, D**, **Z = CN, F** — all
  combinations except **`1D-CN`** and **`1D-F`**.
- Total: **40** `SYSTEM_TYPE_KEYS` (38 multicell + `0A` / `0B`).

---

## Spanwise step-down cascade

Construction narrative:
[`../layout/README.md`](../../layout/README.md).

Web positions are anchored spanwise; webs that no longer fit are **dropped**,
not slid. `effective_layout_for_web_count` / `step_down_cascade` live in
`layout/system_layout.py`.

| X (webs) | Y = A | Y = B | Y = C | Y = D |
|---|---|---|---|---|
| 5 | `5A-CN` `5A-F` | `5B-CN` `5B-F` | `5C-CN` `5C-F` | `5D-CN` `5D-F` |
| 4 | `4A-CN` `4A-F` | `4B-CN` `4B-F` | `4C-CN` `4C-F` | `4D-CN` `4D-F` |
| 3 | `3A-CN` `3A-F` | `3B-CN` `3B-F` | `3C-CN` `3C-F` | `3D-CN` `3D-F` |
| 2 | `2A-CN` `2A-F` | `2B-CN` `2B-F` | `2C-CN` `2C-F` | `2D-CN` `2D-F` |
| 1 | `1A-CN` `1A-F` | `1B-CN` `1B-F` | `1C-CN` `1C-F` | — (no `1D`) |
| 0 | `0A` | `0B` | — | — |

Cascade floors at **X = 1**:

```text
A:  nA-Z → … → 1A-Z
B:  nB-Z → … → 1B-Z
C:  nC-Z → … → 1C-Z
D:  nD-Z → … → 2D-Z → 1C-Z
```

`0A` / `0B` are selectable but not degradation targets. If no web can be
anchored, `InfeasibleSectionGeometry` is raised. **Y = D** hands to **Y = C**
at **X = 1**. Orientation **Z** is carried until **X = 0**.

---

## Reference axes and kwargs (advanced)

These axes are not interchangeable.

| Axis / locus | Role |
|--------------|------|
| **Pitching axis** | Stacking / twist reference; often ~⅓·c from LE. Not elastic axis or shear centre. |
| **Max-thickness locus** | Optional **Y = B** anchor from `AirfoilSDF.thickness_distribution`. |
| **Elastic / centroid / shear centre** | Not used for cap placement in v1. |
| **Web-local reference** | **Y = C** caps centred on each web chord station ([SG-013](../bug_log.md#sg-013)). |

| `structural_family` | Behaviour |
|---------------------|-----------|
| `"D"` *(default)* | Continuous caps between outer webs |
| `"A"` | No spar-cap solids |
| `"B"` | One band; `fixed_cap_anchor` |
| `"C"` | Union of per-web bands |

Kwargs subset: `fixed_cap_anchor`, `pitch_fraction_of_chord_from_le`,
`fixed_cap_chord_half_width`, `discrete_cap_chord_half_width`.

Older docs used **Y = A/B/C** with **C** = box spar; the extended scheme uses
**B**/**C** for discrete patterns and **D** for continuous box.

---

<a id="downstream-handoff"></a>
## Downstream handoff

Gated midlines emit twisted `ShellMidlineStrip` polylines. Prefer
`shell_mesh_inputs_from_geometry_result(result)`. Bind orders zones skin → caps
→ webs; do not assume `panels[i]` aligns with `midlines[i]` when both exist.

---

<a id="caveats-v1"></a>
## Caveats (v1)

- Cap placement is not from shear-centre or elastic-axis coordinates.
- Max-thickness anchor uses sampled `thickness_distribution` only.
- `0A` / `0B` are zero-web hosts: skin ± caps midlines hand off without web junctions.
