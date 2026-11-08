# section examples

Numbered worked examples plus print-first teaching scripts.

↑ [examples](../README.md) · [docs/section](../../docs/section/README.md)

## Numbered

| # | Folder | API |
|---|---|---|
| 01 | [01_layout_catalogue](01_layout_catalogue/) | `SYSTEM_TYPE_KEYS` |
| 02 | [02_section_properties](02_section_properties/) | gated `"properties"` |
| 03 | [03_mid_surface_and_cells](03_mid_surface_and_cells/) | walls, cells, `mid_surface()` |

## Teaching

```bash
python examples/section/teaching/01_first_station.py
```

1. [01_first_station.py](teaching/01_first_station.py) — one twin-web station, all gates
2. [02_read_properties.py](teaching/02_read_properties.py) — `properties_map` table
3. [03_resolve_then_size.py](teaching/03_resolve_then_size.py) — resolve once, two thickness trials
4. [04_spanwise_web_anchors.py](teaching/04_spanwise_web_anchors.py) — anchored webs vs local fractions
5. [05_midline_handoff.py](teaching/05_midline_handoff.py) — print `ShellMidlineStrip` fields (no mesher)
6. [06_consumer_seams.py](teaching/06_consumer_seams.py) — `ThicknessLaw`, named edges, size fields, merged cap+skin

Downstream packages live **outside** this tree. Pin
`implicit-geometry-kernel>=0.8.2` if you read C-family indexed pad keys
(`spar_cap_*_0`) or mapped loft faces; `>=0.8.1` is enough for hybrid C-family
areas and `apply_clip` only (`apply_shapely` is gone). Copy
`teaching/06_consumer_seams.py`. Check a local checkout with
`python scripts/check_consumer.py PATH/TO/CONSUMER`.

Guides: [quickstart](../../docs/section/guides/quickstart.md),
[geometry_to_mesh](../../docs/section/guides/geometry_to_mesh.md).
Gallery: [../demos/section/](../demos/section/).
