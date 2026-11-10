# teaching

Print-first scripts for `SectionGeometry` (`Airfoil.naca4`). Plot-heavy work lives under
[../../demos/section/](../../demos/section/).

↑ [section examples](../README.md)

```bash
python examples/section/teaching/01_first_station.py
```

1. [01_first_station.py](01_first_station.py) — one twin-web station, all gates
2. [02_read_properties.py](02_read_properties.py) — `properties_map` table
3. [03_resolve_then_size.py](03_resolve_then_size.py) — resolve once, two thickness trials
4. [04_spanwise_web_anchors.py](04_spanwise_web_anchors.py) — anchored webs vs local fractions
5. [05_midline_handoff.py](05_midline_handoff.py) — print `ShellMidlineStrip` / `LaminateOffsetRegion` (no mesher)
6. [06_consumer_seams.py](06_consumer_seams.py) — `ThicknessLaw`, named edges, size fields, merged cap+skin

After `05`, feature ids on a cut: [blade/11](../../blade/11_feature_station_network/).
DTO walkthrough: [contract/01](../../contract/01_three_products_dtos/).
