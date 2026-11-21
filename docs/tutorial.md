# Getting started

From two aerofoil stations to a shell mesh in five steps. Every code block on this
page runs in the test suite (`tests/test_tutorial.py`), so it stays correct.

Install the kernel and GEOS first ([Install](https://github.com/samhughes-fea/implicit-geometry-kernel#install)).

## 1. Stations

A station is an aerofoil placed along the span: position `z`, chord, twist and
pitch axis, all in metres and degrees.

```python
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station

stations = [
    Station(0.0, 2.0, Airfoil.naca4("0024"), twist=10.0, pitch_axis=0.3),
    Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
]
```

## 2. A blade

`BladeGeometry` lofts the stations and gives them a structural layout. The code
`2D-F-LC-0R` asks for two shear webs with flat spar caps and no ribs (see the
[catalogue](catalogue/README.md)). Gauges are in metres.

```python
geo = BladeGeometry(
    stations,
    "2D-F-LC-0R",
    skin_thickness=0.003,
    spar_cap_thickness=0.008,
    shear_web_thickness=0.002,
    web_positions=(0.18, 0.48),
    pivot_fraction=0.3,
)
```

## 3. A section

`section_at(z)` cuts the blade at one span position. The cut lists its walls (skin,
caps, webs) as solids with their thickness.

```python
cut = geo.section_at(4.0)
kinds = sorted({wall.kind for wall in cut.walls})
print(kinds)  # ['cap', 'skin', 'web']
```

## 4. The shell graph

The shell graph is what a shell model needs from the section: the skin mid-surface
loop and the web mid-lines, split into partitions that each carry a thickness.

```python
graph = cut.shell_graph()
print(len(graph.loop), "skin runs,", len(graph.webs), "webs")
```

## 5. The shell mesh

`shell_mesh` lofts the shell graphs into a conforming mesh: quads on the skin and
webs, triangles on ribs and the tip closure. Each element carries its surface,
partition and thickness.

```python
from implicit_geometry_kernel.mesh import shell_mesh

mesh = shell_mesh(geo, h=0.25)
print(mesh.nodes.shape[0], "nodes,", len(mesh.elements), "elements")
print(sorted(set(mesh.surface)))  # ['skin', 'tip_closure', 'web']
```

## Next

- [Architecture](architecture.md): how the packages fit together.
- [Worked examples](https://github.com/samhughes-fea/implicit-geometry-kernel/tree/main/examples):
  scripts with figures for every feature.
- [API reference](api/index.md): every public name.
