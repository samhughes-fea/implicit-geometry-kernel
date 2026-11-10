"""Authored WallRegionGraph from catalogue layouts and windIO-shaped tables."""

import numpy as np

from implicit_geometry_kernel import (
    Airfoil,
    Anchor,
    BladeConcept,
    BladeGeometry,
    Station,
    WallRegionGraph,
    transverse_members_from_windio,
)
from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout


def _tapered(key: str, **kwargs):
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, key, 0.003, pivot_fraction=0.3, spar_cap_thickness=0.008,
                         shear_web_thickness=0.002, **kwargs)


def test_family_c_regions_are_per_web_pads_not_a_box():
    graph = WallRegionGraph.from_section_layout("2C-F-LC-0R")
    assert graph.cap_mode == "pads"
    assert len(graph.cap_windows) == 2
    su = [r for r in graph.regions if r.side == "suction" and r.is_cap]
    assert su
    assert ("web_0", "web_1") not in {(r.start, r.end) for r in su}
    fracs = dict(graph.evaluate_anchors(0.0))
    for r in su:
        mid = 0.5 * (fracs[r.start] + fracs[r.end])
        assert any(lo - 1e-9 <= mid <= hi + 1e-9 for lo, hi in graph.cap_windows)


def test_family_d_cap_still_spans_outer_webs():
    graph = WallRegionGraph.from_section_layout("2D-F-LC-0R")
    assert graph.cap_mode == "box"
    su = [r for r in graph.regions if r.side == "suction" and r.is_cap]
    assert [(r.start, r.end) for r in su] == [("web_0", "web_1")]


def test_family_b_cap_is_the_fixed_band():
    graph = WallRegionGraph.from_section_layout("2B-F-LC-0R")
    assert graph.cap_mode == "band"
    su = [r for r in graph.regions if r.side == "suction" and r.is_cap]
    assert any({r.start, r.end} == {"cap_lo", "cap_hi"} for r in su)


def test_overlapping_family_c_pads_merge_into_one_window():
    layout = SectionLayout(
        web_chord_fracs=(0.30, 0.34),
        orientation="flapwise",
        cap_placement=CapPlacementRule.per_web_pads(),
        structural_family="C",
    )
    graph = WallRegionGraph.from_section_layout(layout)
    assert graph.cap_mode == "pads"
    assert len(graph.cap_windows) == 1
    su = [r for r in graph.regions if r.side == "suction" and r.is_cap]
    assert any(r.start.startswith("web_") and r.end.startswith("web_") for r in su)


def test_family_c_one_web_marks_both_pad_halves():
    graph = WallRegionGraph.from_section_layout("1C-F-LC-0R")
    su = [r for r in graph.regions if r.side == "suction" and r.is_cap]
    assert len(su) == 2
    ends = {(r.start, r.end) for r in su}
    assert ("web_0", "cap_hi") in ends
    assert ("cap_lo", "web_0") in ends


def test_from_section_layout_stable_names_catalogue_keys():
    for key in ("0A-LC-0R", "1C-F-LC-0R", "2D-F-LC-0R"):
        graph = WallRegionGraph.from_section_layout(key)
        names = {r.name for r in graph.regions}
        assert "te" in names
        zs = np.linspace(0.0, 10.0, 21)
        sigs = [tuple(r.name for r in graph.regions_at(float(z))) for z in zs]
        assert len(set(sigs)) == 1, (key, set(sigs))
        geo = _tapered(key)
        assert {r.name for r in geo.region_graph.regions} == names
        authored = BladeConcept(key).region_graph()
        assert tuple(r.name for r in authored.regions) == tuple(r.name for r in graph.regions)


def test_te_merge_does_not_spawn_te_adjacent_region():
    layout = SectionLayout.from_catalogue_key("2D-F-LC-0R")
    graph = WallRegionGraph.from_section_layout(layout, te_arc=0.05, min_gap=0.01)
    near = WallRegionGraph(
        graph.anchors + (Anchor("te_ghost", 0.99),),
        graph.regions,
        graph.webs,
        te_arc=0.05,
        min_gap=0.01,
    )
    names = {n for n, _ in near.evaluate_anchors(1.0)}
    assert "te_ghost" not in names
    assert "te" in names


def test_circular_root_region_count_stable():
    stations = [
        Station(0.0, 2.0, Airfoil.ellipse(1.0), pitch_axis=0.5),
        Station(2.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                        shear_web_thickness=0.002, web_positions=(0.18, 0.48), pivot_fraction=0.3)
    counts = [len(geo.region_graph.regions_at(z)) for z in np.linspace(0.0, 10.0, 15)]
    assert len(set(counts)) == 1


def _iea_shaped_structure(*, moving_cap=False, web_drop=False, bulkheads=False):
    g = [0.0, 0.5, 1.0]
    cap_lo = [0.25, 0.28, 0.22] if moving_cap else [0.25, 0.25, 0.25]
    cap_hi = [0.55, 0.52, 0.58] if moving_cap else [0.55, 0.55, 0.55]
    structure = {
        "webs": [
            {
                "name": "fore_web",
                "start_nd_arc": {"grid": g, "values": [0.30, 0.30, 0.30]},
                "end_nd_arc": {"grid": g, "values": [0.70, 0.70, 0.70]},
                **({"start_nd_span": 0.1, "end_nd_span": 0.9} if web_drop else {}),
            },
            {
                "name": "rear_web",
                "start_nd_arc": {"grid": g, "values": [0.20, 0.20, 0.20]},
                "end_nd_arc": {"grid": g, "values": [0.80, 0.80, 0.80]},
            },
        ],
        "layers": [
            {
                "name": "Shell_skin",
                "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.0, 0.0]},
                "end_nd_arc": {"grid": [0.0, 1.0], "values": [1.0, 1.0]},
                "thickness": {"grid": [0.0, 1.0], "values": [0.004, 0.003]},
            },
            {
                "name": "Spar_Cap_SS",
                "start_nd_arc": {"grid": g, "values": cap_lo},
                "end_nd_arc": {"grid": g, "values": cap_hi},
                "thickness": {"grid": [0.0, 1.0], "values": [0.01, 0.008]},
            },
        ],
    }
    if bulkheads:
        structure["bulkheads"] = [
            {"name": "rib_root", "nd_span": 0.1, "thickness": 0.02},
            {"name": "rib_mid", "z": 50.0, "thickness": 0.015},
        ]
    return structure


def test_windio_graph_names_do_not_flicker_with_moving_cap():
    z_of_s = np.column_stack([np.linspace(0, 1, 5), np.linspace(0, 100, 5)])
    graph = WallRegionGraph.from_windio_structure(_iea_shaped_structure(moving_cap=True), z_of_s)
    sigs = [tuple(r.name for r in graph.regions_at(z)) for z in np.linspace(0, 100, 31)]
    assert len(set(sigs)) == 1


def test_transverse_members_tip_closure_and_yaml_bulkheads():
    z_of_s = np.column_stack([[0.0, 1.0], [0.0, 80.0]])
    members = transverse_members_from_windio(_iea_shaped_structure(), z_of_s)
    assert [m.name for m in members] == ["tip_closure"]
    assert members[0].z == 80.0
    extra = transverse_members_from_windio(_iea_shaped_structure(bulkheads=True), z_of_s)
    names = {m.name for m in extra}
    assert names == {"tip_closure", "rib_root", "rib_mid"}
    assert {round(m.z, 6) for m in extra if m.name != "tip_closure"} == {8.0, 50.0}
