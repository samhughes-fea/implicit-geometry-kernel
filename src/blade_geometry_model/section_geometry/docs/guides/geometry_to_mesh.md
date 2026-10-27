# Geometry to mesh

Turn gated midlines into `section_mesh` inputs. Overview:
[`../README.md`](../../README.md). Runnable twin:
[`../examples/05_handoff_to_mesh.py`](../../examples/05_handoff_to_mesh.py).

↑ [`section_geometry`](../../README.md) · [`docs index`](../README.md)

## Table of contents

1. [Pipeline](#pipeline)
2. [Required geometry gates](#required-geometry-gates)
3. [Adapter](#adapter)
4. [What `ShellMeshInputs` carries](#what-shellmeshinputs-carries)
5. [Optional coarse strip mesh](#optional-coarse-strip-mesh)
6. [What next](#what-next)

---

## Pipeline

```text
section_geometry (midlines) → shell_mesh_inputs_from_geometry_result
                           → section_mesh (strips)
                           → section_bind (materials; after mesh)
```

This package stays material-bind agnostic: handoff is wall centrelines and
thicknesses only.

---

## Required geometry gates

Include `"midlines"` (or `"cell_wall_loops"`, which also fills midlines):

```python
result = SectionGeometryRunner(
    layout=resolve_system_type("2D-F"),
    airfoil=airfoil,
    skin_thickness=0.004,
    spar_cap_thickness=0.012,
    shear_web_thickness=0.005,
    outputs=("midlines",),
).execute().get_results()
```

Layouts `"0A"` / `"0B"` are valid shell-mesh stations (skin ± caps; no web
junctions). See
[`ADR_ZERO_WEB_SHELL_MESH.md`](../reference/ADR_ZERO_WEB_SHELL_MESH.md) and
[`feasibility_and_drops.md`](feasibility_and_drops.md#zero-web-keys).

---

## Adapter

```python
from blade_structure_utilities.mesh_utilities.section_mesh import (
    shell_mesh_inputs_from_geometry_result,
)

shell = shell_mesh_inputs_from_geometry_result(result, layout_key="2D-F")
```

Raises `ValueError` if `midline_strips` (or chord / twist / LE/TE) are missing.

---

## What `ShellMeshInputs` carries

| Field | Role |
|-------|------|
| `midlines` | `ShellMidlineStrip` list (skin / web / cap polylines + `thickness_m`) |
| `chord_m`, `twist_rad` | Station scalars |
| `leading_edge_…`, `trailing_edge_…` | Twisted LE/TE markers in the blade frame |
| `layout_key` | Optional catalogue key for downstream logging |

Strip `kind` / `label` match the anatomy table in
[`anatomy_on_a_figure.md`](anatomy_on_a_figure.md)
(same folder as this guide).

---

## Optional coarse strip mesh

```python
from blade_structure_utilities.mesh_utilities.section_mesh import (
    build_section_strip_mesh,
)

mesh = build_section_strip_mesh(shell, n_elements_per_panel=4, endpoint_tol=1e-6)
```

Example flag: `python .../05_handoff_to_mesh.py --mesh`.

Figure demos for strip meshes live under
[`../../../mesh_utilities/section_mesh/visualisation/demos/`](../../../../mesh_utilities/section_mesh/visualisation/demos).

---

## What next

| Goal | Go to |
|------|--------|
| Bind materials on the mesh | [`section_bind`](../../../../material_utilities/section_bind/README.md) · [`examples/`](../../../../material_utilities/section_bind/examples/README.md) |
| Shape `SectionForceVector` | [`section_loads`](../../../../load_utilities/section_loads/README.md) · [`examples/`](../../../../load_utilities/section_loads/examples/README.md) · [`quickstart`](../../../../load_utilities/section_loads/docs/quickstart.md) |
| Spine ownership | [`COEXISTENCE_IMPLICIT_SHELL`](../../../../docs/COEXISTENCE_IMPLICIT_SHELL.md) |
| Midline emit details | [`../midlines/`](../../midlines/README.md) |
