# 01 · Transverse rib face on the shell graph

A `TransverseMember` is a rib, bulkhead or frame at one span station. It is not a
section wall type and does not change the catalogue letter. `BladeGeometry` stores
the members; `section_at` at that `z` attaches closed rib faces on
`shell_graph().faces` from the cell-wall loops. Other stations have no faces.
`forced_stations()` lists the member `z` values as a mesh hint.

```bash
python examples/structure/01_transverse_rib/transverse_rib.py
```

```python
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, TransverseMember
# A rib is a TransverseMember on BladeGeometry, not a BladeConcept catalogue letter.

stations = [Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
            Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3)]
rib = TransverseMember("rib_mid", z=5.0, thickness=0.012)
blade = BladeGeometry(stations, "2D-F-LC-0R", 0.003, pivot_fraction=0.3, spar_cap_thickness=0.008,
                      shear_web_thickness=0.002, web_positions=(0.18, 0.48),
                      transverse=(rib,))
faces = blade.section_at(5.0).shell_graph().faces
assert blade.section_at(2.0).shell_graph().faces == ()
```
