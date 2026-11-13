"""Geometric integrals include a spar cap on its own mid-line."""

import numpy as np

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.cut import SectionWall
from implicit_geometry_kernel.blade.geometric import geometric_integrals


def _wall(xy, thickness, kind="skin"):
    pts = np.asarray(xy, dtype=float)
    pts = np.column_stack([pts, np.zeros(len(pts))])
    return SectionWall(label=kind, kind=kind, points=pts, thickness=thickness)


def test_a_strip_offset_from_the_axis_adds_area_times_offset_squared():
    on_axis = geometric_integrals([_wall([[0.0, 0.0], [1.0, 0.0]], 0.02)])
    offset = geometric_integrals([_wall([[0.0, 0.05], [1.0, 0.05]], 0.02)])
    assert on_axis.area == offset.area
    assert abs(on_axis.ixx_origin) < 1e-15
    assert abs(offset.ixx_origin - offset.area * 0.05 ** 2) < 1e-15
    assert abs(offset.centroid[1] - 0.05) < 1e-15


def test_a_station_keeps_the_cap_offset_in_the_moments():
    foil = Airfoil.naca4("0024")
    stations = [
        Station(0.0, 1.0, foil, pitch_axis=0.3),
        Station(4.0, 1.0, foil, pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.004,
        web_positions=(0.2, 0.55),
        pivot_fraction=0.3,
    )
    cut = geo.section_at(2.0)
    full = geometric_integrals(cut)
    caps = [w for w in cut.walls if w.kind == "cap"]
    assert caps
    without_caps = geometric_integrals([w for w in cut.walls if w.kind != "cap"])
    assert full.area > without_caps.area
    collapsed = []
    for wall in cut.walls:
        if wall.kind != "cap":
            collapsed.append(wall)
            continue
        pts = np.array(wall.points, dtype=float)
        pts[:, 1] = 0.0
        collapsed.append(SectionWall(wall.label, wall.kind, pts, wall.thickness))
    flat = geometric_integrals(collapsed)
    assert abs(full.ixx_origin - flat.ixx_origin) > 1e-8
