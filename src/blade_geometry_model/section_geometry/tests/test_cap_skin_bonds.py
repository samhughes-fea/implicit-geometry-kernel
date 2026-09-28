"""Caps are bonded to their skin by an emitted vertex correspondence (SG-017)."""

from __future__ import annotations

import itertools
import warnings

import numpy as np
import pytest

LAYOUTS = ["1C-F", "1C-CN", "2D-F", "2D-CN", "3D-F", "3C-F", "4D-F", "4D-CN", "2B-F"]
CASES = list(itertools.product(LAYOUTS, (0.0, 0.6, 1.2), (1.0, 0.2)))


def _build(key: str, twist_rad: float, chord: float, **kwargs):
    from blade_structure_utilities.load_utilities.section_loads.visualisation.demos._common import (
        build_midlines_for_key,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_midlines_for_key(
            key,
            chord=chord,
            twist_rad=twist_rad,
            skin_thickness=0.004 * chord,
            spar_cap_thickness=0.01 * chord,
            shear_web_thickness=0.004 * chord,
            **kwargs,
        )


def _on_polyline(point: np.ndarray, polyline: np.ndarray) -> float:
    a, b = polyline[:-1], polyline[1:]
    ab = b - a
    t = np.clip(np.einsum("ij,ij->i", point - a, ab) / np.einsum("ij,ij->i", ab, ab), 0.0, 1.0)
    return float(np.min(np.hypot(*(a + t[:, None] * ab - point).T)))


@pytest.mark.parametrize("key, twist_rad, chord", CASES)
def test_every_cap_is_bonded_one_to_one_along_its_whole_length(key, twist_rad, chord):
    strips, _loops, _ = _build(key, twist_rad, chord)
    by_label = {s.label: s for s in strips}
    caps = [s for s in strips if s.kind == "cap"]
    assert caps
    for cap in caps:
        assert cap.bonds, f"{cap.label} has no bond"
        cap_points = np.asarray(cap.midline_twisted_pitch_axis_blade_frame)
        covered: set[int] = set()
        for bond in cap.bonds:
            host = np.asarray(by_label[bond.host_label].midline_twisted_pitch_axis_blade_frame)
            assert by_label[bond.host_label].kind == "skin"
            cap_index = [int(np.argmin(np.hypot(*(cap_points - p).T))) for p in bond.points]
            host_index = [int(np.argmin(np.hypot(*(host - p).T))) for p in bond.host_points]
            assert np.allclose(cap_points[cap_index], bond.points, atol=1e-12, rtol=0)
            assert np.allclose(host[host_index], bond.host_points, atol=1e-12, rtol=0)
            assert set(np.abs(np.diff(cap_index))) == {1}, "bonded cap vertices must be consecutive"
            assert set(np.abs(np.diff(host_index))) == {1}, "bonded skin vertices must be consecutive"
            for a, b in zip(cap_index[:-1], cap_index[1:]):
                segment = min(a, b)
                assert segment not in covered, "a cap segment is bonded twice"
                covered.add(segment)
        assert covered == set(range(len(cap_points) - 1)), f"{cap.label} is only partly bonded"


@pytest.mark.parametrize("key", ["2D-F", "3C-F", "4D-F"])
def test_bonding_moves_no_geometry(key):
    """Every pre-bond vertex survives and every added vertex lies on its original polyline."""
    from blade_structure_utilities.geometry_utilities.section_geometry.topology import shell_strip_topology

    captured = {}
    original = shell_strip_topology.bond_caps_to_skin

    def spy(strips):
        captured["before"] = {s.label: np.asarray(s.midline_twisted_pitch_axis_blade_frame) for s in strips}
        result = original(strips)
        captured["after"] = {s.label: np.asarray(s.midline_twisted_pitch_axis_blade_frame) for s in result}
        return result

    shell_strip_topology.bond_caps_to_skin = spy
    try:
        _build(key, 0.6, 0.2)
    finally:
        shell_strip_topology.bond_caps_to_skin = original
    for label, before in captured["before"].items():
        after = captured["after"][label]
        for point in before:
            assert np.min(np.hypot(*(after - point).T)) == 0.0, f"{label}: a vertex moved"
        for point in after:
            assert _on_polyline(point, before) < 1e-12, f"{label}: an added vertex is off the polyline"


@pytest.mark.parametrize("key, twist_rad, chord", CASES)
def test_the_thin_wall_kernel_is_exact_on_bonded_worker_sections(key, twist_rad, chord):
    """With the emitted bonds, every unit load integrates back exactly (TW-001)."""
    from blade_structure_utilities.geometry_utilities.section_thin_wall_model import (
        ThinWall,
        WallBond,
        axial_field,
        axial_stiffness,
        build_wall_graph,
        contour_samples,
        integrate_samples,
        shear_flow,
    )

    strips, loops, _ = _build(key, twist_rad, chord)
    walls = [
        ThinWall(
            s.label,
            s.kind,
            s.midline_twisted_pitch_axis_blade_frame,
            7.0e10 * s.thickness_m,
            0.0 if s.kind == "cap" else 2.7e10 * s.thickness_m,
            closed=bool(s.closed),
            carries_shear=s.kind != "cap",
        )
        for s in strips
    ]
    bonds = [WallBond(s.label, b.host_label, b.points, b.host_points) for s in strips for b in s.bonds]
    graph = build_wall_graph(walls, [loop.wall_labels for loop in loops], bonds=bonds)
    reference = (0.01 * chord, -0.02 * chord)
    stiffness = axial_stiffness(graph, reference)
    names = ("N", "Vy", "Vz", "My", "Mz", "T")
    for row, name in enumerate(names):
        loads = dict.fromkeys(names, 0.0)
        loads[name] = 1.0
        axial = axial_field(stiffness, N=loads["N"], My=loads["My"], Mz=loads["Mz"])
        flow = shear_flow(graph, stiffness, Vz=loads["Vz"], Vy=loads["Vy"], T=loads["T"])
        got = integrate_samples(contour_samples(graph, axial, flow), reference).as_array()
        assert np.max(np.abs(got - np.eye(6)[row])) < 1e-9, (name, got)
