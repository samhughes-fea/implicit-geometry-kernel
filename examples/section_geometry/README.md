# examples

Short, runnable teaching scripts for `section_geometry`. Plot-heavy figure
generation lives under [`../visualisation/demos/`](../visualisation/demos/README.md).

↑ [`section_geometry`](../README.md) · [`docs index`](../docs/README.md)

## Table of contents

1. [Run from repo root](#run-from-repo-root)
2. [Newcomer path](#newcomer-path)
3. [Scripts](#scripts)
4. [See also](#see-also)

---

## Run from repo root

```bash
python blade_structure_utilities/geometry_utilities/section_geometry/examples/<script>.py
```

Or install editable (`pip install -e .`) and run the same path.

## Newcomer path

1. [`01_first_station.py`](01_first_station.py) — one twin-web station, all gates
2. [`02_read_properties.py`](02_read_properties.py) — `properties_map` table
3. [`03_resolve_then_size.py`](03_resolve_then_size.py) — resolve once, two thickness trials
4. [`04_spanwise_web_anchors.py`](04_spanwise_web_anchors.py) — anchored webs vs local fractions
5. [`05_handoff_to_mesh.py`](05_handoff_to_mesh.py) — midlines → `ShellMeshInputs`

## Scripts

| Script | Intent |
|--------|--------|
| [`01_first_station.py`](01_first_station.py) | NACA 2412 + `"2D-F"`; print labels / strips / loops (`--plot` optional) |
| [`02_read_properties.py`](02_read_properties.py) | Area / centroid table from `"properties"` |
| [`03_resolve_then_size.py`](03_resolve_then_size.py) | `resolve()` once, rebuild with `resolution=` |
| [`04_spanwise_web_anchors.py`](04_spanwise_web_anchors.py) | Tip station: local fractions vs root-anchored drops |
| [`05_handoff_to_mesh.py`](05_handoff_to_mesh.py) | `shell_mesh_inputs_from_geometry_result` (`--mesh` optional) |

## See also

- Progressive guide: [`../docs/guides/quickstart.md`](../docs/guides/quickstart.md)
- Call recipes: [`../docs/guides/caller_patterns.md`](../docs/guides/caller_patterns.md)
- Anatomy on a figure: [`../docs/guides/anatomy_on_a_figure.md`](../docs/guides/anatomy_on_a_figure.md)
- Feasibility / drops: [`../docs/guides/feasibility_and_drops.md`](../docs/guides/feasibility_and_drops.md)
- Geometry → mesh: [`../docs/guides/geometry_to_mesh.md`](../docs/guides/geometry_to_mesh.md)
- Figure demos: [`../visualisation/demos/`](../visualisation/demos/README.md)
