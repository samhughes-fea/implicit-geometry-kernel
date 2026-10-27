# Structural systems gallery

All 40 `SYSTEM_TYPE_KEYS`. Pass any key to `resolve_system_type(...)` /
`SectionGeometryRunner(layout=...)`. Grouped by web count **X** (≠ spar-cap
family **Y**). Decode: [`system_type_xyz_taxonomy.md`](system_type_xyz_taxonomy.md).

↑ [`section_geometry`](../../README.md) · [`../guides/quickstart.md`](../guides/quickstart.md)

Figures: NACA 2412 + demo thicknesses from `structural_family_grids_demo.py`
([`../visualisation/output/structural_family_grids/`](../../visualisation/output/structural_family_grids);
gitignored — regenerate to preview; see
[`../visualisation/demos/ARTIFACTS.md`](../../visualisation/demos/ARTIFACTS.md)).
Production builds set `skin_thickness` /
`spar_cap_thickness` / `shear_web_thickness` (optional `web_x_positions_m`).

## Table of contents

1. [Web count X = 0](#web-x-0)
2. [Web count X = 1](#web-x-1)
3. [Web count X = 2](#web-x-2)
4. [Web count X = 3](#web-x-3)
5. [Web count X = 4](#web-x-4)
6. [Web count X = 5](#web-x-5)
7. [See also](#see-also)

**Regenerate** (from repo root; full artefact list in
[`ARTIFACTS.md`](../../visualisation/demos/ARTIFACTS.md)):

```bash
python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/structural_family_grids_demo.py
```

---

<a id="web-x-0"></a>
## Web count X = 0 — no shear webs

Skin (optional single-band caps). `0A` / `0B` remain catalogue keys; shell-mesh
handoff is skin ± caps with no web junctions
([`ADR_ZERO_WEB_SHELL_MESH.md`](ADR_ZERO_WEB_SHELL_MESH.md)).

| Key | Meaning |
|-----|---------|
| `0A` | Outer skin only (no spar caps) |
| `0B` | Skin + one upper/lower spar-cap band |

### `0A`

![0A](../../visualisation/output/structural_family_grids/0A.png)

### `0B`

![0B](../../visualisation/output/structural_family_grids/0B.png)

---

<a id="web-x-1"></a>
## Web count X = 1 — one shear web

Single web plus spar-cap family **Y** = **A** / **B** / **C**. Continuous box
(**Y** = **D**) is not registered for one web.

| Key | Meaning |
|-----|---------|
| `1A-CN` / `1A-F` | One web, no spar caps; CN or flapwise |
| `1B-CN` / `1B-F` | One web + single fixed cap band |
| `1C-CN` / `1C-F` | One web + discrete cap band at the web |

### `1A-CN`

![1A-CN](../../visualisation/output/structural_family_grids/1A-CN.png)

### `1A-F`

![1A-F](../../visualisation/output/structural_family_grids/1A-F.png)

### `1B-CN`

![1B-CN](../../visualisation/output/structural_family_grids/1B-CN.png)

### `1B-F`

![1B-F](../../visualisation/output/structural_family_grids/1B-F.png)

### `1C-CN`

![1C-CN](../../visualisation/output/structural_family_grids/1C-CN.png)

### `1C-F`

![1C-F](../../visualisation/output/structural_family_grids/1C-F.png)

---

<a id="web-x-2"></a>
## Web count X = 2 — two shear webs

Typical twin-web / torsion-box station. Includes continuous box (**Y** = **D**).

| Key | Meaning |
|-----|---------|
| `2A-*` | Two webs, no spar caps |
| `2B-*` | Two webs + one fixed cap band |
| `2C-*` | Two webs + discrete caps per web |
| `2D-*` | Two webs + continuous box caps between outer webs |

### `2A-CN`

![2A-CN](../../visualisation/output/structural_family_grids/2A-CN.png)

### `2A-F`

![2A-F](../../visualisation/output/structural_family_grids/2A-F.png)

### `2B-CN`

![2B-CN](../../visualisation/output/structural_family_grids/2B-CN.png)

### `2B-F`

![2B-F](../../visualisation/output/structural_family_grids/2B-F.png)

### `2C-CN`

![2C-CN](../../visualisation/output/structural_family_grids/2C-CN.png)

### `2C-F`

![2C-F](../../visualisation/output/structural_family_grids/2C-F.png)

### `2D-CN`

![2D-CN](../../visualisation/output/structural_family_grids/2D-CN.png)

### `2D-F`

![2D-F](../../visualisation/output/structural_family_grids/2D-F.png)

---

<a id="web-x-3"></a>
## Web count X = 3 — three shear webs

Same Y/Z grid as X = 2, with a third web (more cells along chord).

### `3A-CN`

![3A-CN](../../visualisation/output/structural_family_grids/3A-CN.png)

### `3A-F`

![3A-F](../../visualisation/output/structural_family_grids/3A-F.png)

### `3B-CN`

![3B-CN](../../visualisation/output/structural_family_grids/3B-CN.png)

### `3B-F`

![3B-F](../../visualisation/output/structural_family_grids/3B-F.png)

### `3C-CN`

![3C-CN](../../visualisation/output/structural_family_grids/3C-CN.png)

### `3C-F`

![3C-F](../../visualisation/output/structural_family_grids/3C-F.png)

### `3D-CN`

![3D-CN](../../visualisation/output/structural_family_grids/3D-CN.png)

### `3D-F`

![3D-F](../../visualisation/output/structural_family_grids/3D-F.png)

---

<a id="web-x-4"></a>
## Web count X = 4 — four shear webs

Higher web count; default chord fractions come from the layout registry.

### `4A-CN`

![4A-CN](../../visualisation/output/structural_family_grids/4A-CN.png)

### `4A-F`

![4A-F](../../visualisation/output/structural_family_grids/4A-F.png)

### `4B-CN`

![4B-CN](../../visualisation/output/structural_family_grids/4B-CN.png)

### `4B-F`

![4B-F](../../visualisation/output/structural_family_grids/4B-F.png)

### `4C-CN`

![4C-CN](../../visualisation/output/structural_family_grids/4C-CN.png)

### `4C-F`

![4C-F](../../visualisation/output/structural_family_grids/4C-F.png)

### `4D-CN`

![4D-CN](../../visualisation/output/structural_family_grids/4D-CN.png)

### `4D-F`

![4D-F](../../visualisation/output/structural_family_grids/4D-F.png)

---

<a id="web-x-5"></a>
## Web count X = 5 — five shear webs

Maximum registered web count in the layout catalogue
([`../layout/system_layout.py`](../../layout/system_layout.py)).

### `5A-CN`

![5A-CN](../../visualisation/output/structural_family_grids/5A-CN.png)

### `5A-F`

![5A-F](../../visualisation/output/structural_family_grids/5A-F.png)

### `5B-CN`

![5B-CN](../../visualisation/output/structural_family_grids/5B-CN.png)

### `5B-F`

![5B-F](../../visualisation/output/structural_family_grids/5B-F.png)

### `5C-CN`

![5C-CN](../../visualisation/output/structural_family_grids/5C-CN.png)

### `5C-F`

![5C-F](../../visualisation/output/structural_family_grids/5C-F.png)

### `5D-CN`

![5D-CN](../../visualisation/output/structural_family_grids/5D-CN.png)

### `5D-F`

![5D-F](../../visualisation/output/structural_family_grids/5D-F.png)

---

## See also

- [`../guides/quickstart.md`](../guides/quickstart.md) — progressive call guide
- [`system_type_xyz_taxonomy.md`](system_type_xyz_taxonomy.md) — id decode, cascade
- [`../guides/caller_patterns.md`](../guides/caller_patterns.md) — sizing loop and spanwise webs
- [`../layout/README.md`](../../layout/README.md) — layout resolve and web drop
- Demo script: `visualisation/demos/structural_family_grids_demo.py`
