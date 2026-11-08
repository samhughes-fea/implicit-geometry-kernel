"""Step 6 seams: thickness law, named edges, size fields, merged cap+skin."""

from __future__ import annotations

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.stations import PlyDrop, TransverseMember
from implicit_geometry_kernel.structure.laws import ConstantLaw, ThicknessLaw


def _geo(**kwargs):
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations,
        "2D-F-LC-0R",
        0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48),
        pivot_fraction=0.3,
        **kwargs,
    )


def test_thickness_law_from_span_ignores_s():
    law = ThicknessLaw.from_span(ConstantLaw(0.004))
    t = law(np.array([1.0, 2.0]), np.array([0.1, 0.9]))
    np.testing.assert_allclose(t, 0.004)


def test_blade_thickness_law_matches_thickness_at():
    drop = PlyDrop("d1", 0.2, 0.4, 1.0, 5.0, 0.002, taper=0.0)
    geo = _geo(ply_drops=(drop,))
    law = geo.thickness_law()
    np.testing.assert_allclose(law(3.0, 0.3), geo.thickness_at(3.0, 0.3))
    assert float(law(3.0, 0.3)) == pytest.approx(0.005)


def test_named_edges_include_partitions_and_ply_drops():
    drop = PlyDrop("d1", 0.2, 0.4, 1.0, 5.0, 0.002)
    geo = _geo(ply_drops=(drop,))
    edges = geo.named_edges(4.0)
    names = {(e.name_a, e.name_b) for e in edges}
    assert any("d1:s0" in pair for pair in names)
    assert any("d1:s1" in pair for pair in names)
    assert geo.section_at(4.0).shell_graph().partition_edges()


def test_forced_stations_include_ply_drop_ends():
    drop = PlyDrop("d1", 0.2, 0.4, 1.0, 5.0, 0.002)
    rib = TransverseMember("rib", 3.0, 0.01)
    geo = _geo(ply_drops=(drop,), transverse=(rib,))
    assert geo.forced_stations() == (1.0, 3.0, 5.0)


def test_size_fields_keys():
    geo = _geo()
    cut = geo.section_at(4.0)
    p = cut.shell_graph().loop[0].points[len(cut.shell_graph().loop[0].points) // 2]
    out = geo.size_fields(p.reshape(1, 3))
    assert set(out) >= {"thickness", "distance", "curvature", "member"}
    assert out["curvature"].shape == (1,)


def test_merged_cap_and_skin_uses_graph_runs():
    stacked = _geo()
    merged = _geo(merge_cap_and_skin=True)
    a = stacked.section_at(4.0)
    b = merged.section_at(4.0)
    assert len(a.walls) == 10
    assert len(b.walls) == len(a.shell_graph().runs())
    assert any(w.kind == "cap" for w in b.walls)
    cap = next(w for w in b.walls if w.kind == "cap")
    assert cap.thickness == pytest.approx(0.011)
