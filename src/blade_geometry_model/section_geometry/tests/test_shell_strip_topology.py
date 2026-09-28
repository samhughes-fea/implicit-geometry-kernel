"""Tests for geometric shell strip junction topology (section_geometry)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_REPO = Path(__file__).resolve().parents[3]
for _p in (str(_REPO),):
    if _p not in sys.path:
        sys.path.insert(0, _p)


@pytest.fixture(scope="module")
def airfoil():
    from blade_structure_utilities.geometry_utilities.section_geometry.sdf.airfoil import AirfoilSDF

    return AirfoilSDF.from_naca("0012", n_points=100, chord=1.0)


def _strips(airfoil, layout_key: str):
    from blade_structure_utilities.geometry_utilities.section_geometry.api.factory import build_section_view, resolve_system_type
    from blade_structure_utilities.geometry_utilities.section_geometry.midlines.shell_midline_export import build_shell_midline_strips

    layout = resolve_system_type(layout_key)
    section = build_section_view(airfoil, layout, twist_angle_rad=0.0)
    return build_shell_midline_strips(
        section, twist_rad=0.0, n_web_samples=20, n_cap_samples=40
    )


@pytest.mark.parametrize("layout_key", ["2D-F", "2D-CN"])
def test_prepare_shell_midline_junctions_preserves_cap_endpoints(layout_key, airfoil):
    from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation.shell_midline_junction_preparation import (
        ShellMidlineJunctionPreparationOptions,
        prepare_shell_midline_junctions,
        split_skin_at_junctions,
    )

    strips = _strips(airfoil, layout_key)
    split = split_skin_at_junctions(strips)
    ends_split = {
        m.label: (np.asarray(m.midline_twisted_pitch_axis_blade_frame)[0].copy(), np.asarray(m.midline_twisted_pitch_axis_blade_frame)[-1].copy())
        for m in split
        if m.kind == "cap"
    }
    out = prepare_shell_midline_junctions(strips, ShellMidlineJunctionPreparationOptions())
    for m in out:
        if m.kind != "cap":
            continue
        arr = np.asarray(m.midline_twisted_pitch_axis_blade_frame, dtype=float)
        e0, e1 = ends_split[m.label]
        assert np.allclose(arr[0], e0, atol=1e-9)
        assert np.allclose(arr[-1], e1, atol=1e-9)


@pytest.mark.parametrize("layout_key", ["2C-CN", "3D-F"])
def test_topology_inserts_cap_web_junctions(layout_key, airfoil):
    from blade_structure_utilities.geometry_utilities.section_geometry.api.factory import build_section_view, resolve_system_type
    from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation.shell_midline_junction_preparation import (
        ShellMidlineJunctionPreparationOptions,
        prepare_shell_midline_junctions,
    )
    from blade_structure_utilities.mesh_utilities.section_mesh.handoff._section_to_shell_inputs import (
        shell_mesh_inputs_from_section,
    )

    layout = resolve_system_type(layout_key)
    section = build_section_view(airfoil, layout, twist_angle_rad=0.0)
    strips = _strips(airfoil, layout_key)
    topo = prepare_shell_midline_junctions(strips, ShellMidlineJunctionPreparationOptions())
    shell = shell_mesh_inputs_from_section(
        section,
        twist_rad=0.0,
        layout_key=layout_key,
        n_web_samples=20,
        n_cap_samples=80,
    )
    topo_caps = {m.label: np.asarray(m.midline_twisted_pitch_axis_blade_frame).shape[0] for m in topo if m.kind == "cap"}
    shell_caps = {m.label: np.asarray(m.midline_twisted_pitch_axis_blade_frame).shape[0] for m in shell.midlines if m.kind == "cap"}
    assert topo_caps == shell_caps


@pytest.mark.parametrize("layout_key, expect_kinds", [("0A", {"skin"}), ("0B", {"skin", "cap"})])
def test_zero_web_topology_skips_web_junctions(layout_key, expect_kinds, airfoil):
    from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation.shell_midline_junction_preparation import (
        ShellMidlineJunctionPreparationOptions,
        prepare_shell_midline_junctions,
    )

    strips = _strips(airfoil, layout_key)
    assert not any(s.kind == "web" for s in strips)
    topo = prepare_shell_midline_junctions(strips, ShellMidlineJunctionPreparationOptions())
    kinds = {s.kind for s in topo}
    assert kinds == expect_kinds
    assert not any(s.kind == "web" for s in topo)
    for s in topo:
        if s.kind == "skin":
            assert not str(s.label).endswith(":leading_edge")
            assert not str(s.label).endswith(":trailing_edge")


def test_zero_web_cell_loops_single_cavity(airfoil):
    from blade_structure_utilities.geometry_utilities.section_geometry.cell_wall_loops import (
        build_cell_wall_loops,
    )
    from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation.shell_midline_junction_preparation import (
        ShellMidlineJunctionPreparationOptions,
        prepare_shell_midline_junctions,
    )

    strips = _strips(airfoil, "0A")
    topo = prepare_shell_midline_junctions(strips, ShellMidlineJunctionPreparationOptions())
    loops = build_cell_wall_loops(topo)
    assert len(loops) == 1
    assert loops[0].web_labels == ()