"""Public geometry worker facade: gates, skin-only path, prep key."""

from __future__ import annotations

import pytest

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    build_section_geometry,
    build_section_with_thicknesses,
    resolve_system_type,
    section_geometry_prep_key,
)
from blade_structure_utilities.geometry_utilities.section_geometry.api.worker import (
    resolve_station_layout,
)


def test_build_section_with_thicknesses_multicell():
    af = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    layout = resolve_system_type("1B-CN")
    section = build_section_with_thicknesses(
        af, layout, twist_rad=0.0, skin_thickness=0.005, spar_cap_thickness=0.008, shear_web_thickness=0.006
    )
    assert "skin_suction_surface" in section.labels
    assert "skin_pressure_surface" in section.labels
    assert any(
        str(getattr(lbl, "value", lbl)).startswith(("shear_web", "web"))
        for lbl in section.labels
    )


def test_build_section_with_thicknesses_skin_only_public():
    from blade_structure_utilities.geometry_utilities.section_geometry.sections import (
        MultiCellSection,
    )

    af = AirfoilSDF.from_naca("0012", n_points=60, chord=1.0)
    layout = resolve_system_type("0A")
    section = build_section_with_thicknesses(
        af, layout, twist_rad=0.0, skin_thickness=0.004, spar_cap_thickness=0.0, shear_web_thickness=0.0
    )
    assert isinstance(section, MultiCellSection)
    assert section.n_webs == 0
    assert "skin_suction_surface" in section.labels
    assert "skin_pressure_surface" in section.labels
    assert "spar_cap_suction_surface" not in section.labels


def test_build_section_with_thicknesses_0B_is_multicell():
    from blade_structure_utilities.geometry_utilities.section_geometry.sections import (
        MultiCellSection,
    )
    from blade_structure_utilities.geometry_utilities.section_geometry.sections.from_layout import (
        fixed_cap_anchor_from_layout,
    )

    af = AirfoilSDF.from_naca("0012", n_points=80, chord=1.0)
    layout = resolve_system_type("0B")
    assert fixed_cap_anchor_from_layout(layout) == "max_thickness"
    section = build_section_with_thicknesses(
        af, layout, twist_rad=0.0, skin_thickness=0.004, spar_cap_thickness=0.008, shear_web_thickness=0.0
    )
    assert isinstance(section, MultiCellSection)
    assert section.n_webs == 0
    assert "spar_cap_suction_surface" in section.labels
    assert "spar_cap_pressure_surface" in section.labels

def test_build_section_geometry_midlines_gate_skips_props():
    af = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    layout = resolve_system_type("1B-CN")
    result = build_section_geometry(
        layout,
        af,
        skin_thickness=0.005,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.006,
        outputs=("midlines",),
    )
    assert result.properties_map is None
    assert result.midline_strips is not None
    assert len(result.midline_strips) >= 1
    assert result.cell_wall_loops is None
    assert result.chord_m == pytest.approx(1.0)
    assert result.twist_rad == pytest.approx(0.0)
    assert result.leading_edge_twisted_pitch_axis_blade_frame is not None
    assert result.trailing_edge_twisted_pitch_axis_blade_frame is not None
    key = section_geometry_prep_key(
        layout, af, skin_thickness=0.005, spar_cap_thickness=0.008, shear_web_thickness=0.006, outputs=("midlines",)
    )
    assert result.prep_key == key


def test_stage_timings_attached_when_profiling_enabled(monkeypatch):
    monkeypatch.setenv("SECTION_GEOMETRY_PROFILE", "1")
    af = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    layout = resolve_system_type("1B-CN")
    result = build_section_geometry(
        layout,
        af,
        skin_thickness=0.005,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.006,
        outputs=("properties", "midlines"),
    )
    assert result.stage_timings is not None
    assert "assemble" in result.stage_timings
    assert "assemble_skin" in result.stage_timings
    assert "assemble_anchors" in result.stage_timings
    assert "assemble_webs" in result.stage_timings
    assert "assemble_caps" in result.stage_timings
    assert "assemble_cores" in result.stage_timings
    assert "properties" in result.stage_timings
    assert "properties_constructive" in result.stage_timings
    assert "properties_residual" in result.stage_timings
    assert "emit_wall_solids" in result.stage_timings
    assert "emit_unions" in result.stage_timings
    assert "emit_cavity" in result.stage_timings
    assert "emit_cores" in result.stage_timings
    assert "properties_shoelace" in result.stage_timings
    assert "midlines" in result.stage_timings
    assert "topology" in result.stage_timings
    assert result.stage_timings["assemble"] >= 0.0
    assert result.stage_timings["properties"] >= 0.0
    assert result.stage_timings["properties_constructive"] >= 0.0
    assert result.stage_timings["properties_residual"] >= 0.0
    assert result.stage_timings["emit_wall_solids"] >= 0.0
    assert result.stage_timings["properties_shoelace"] >= 0.0


def test_resolve_station_layout_anchors_skip_second_resolve_web_anchors(monkeypatch):
    """G8: layout resolution already ran resolve_web_anchors; build must not redo it."""
    import blade_structure_utilities.geometry_utilities.section_geometry.sections.multicell as mc

    af = AirfoilSDF.from_naca("0012", n_points=120, chord=1.0)
    layout = resolve_system_type("2D-CN")
    anchored = [float(f) * 1.0 for f in layout.web_chord_fracs]
    station = resolve_station_layout(
        layout, af, skin_thickness=0.005, web_x_positions_m=anchored
    )
    assert station.web_y_coords
    assert len(station.web_y_coords) == len(station.web_x_positions_m)

    calls = {"n": 0}
    real = mc.resolve_web_anchors

    def _counting(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    monkeypatch.setattr(mc, "resolve_web_anchors", _counting)
    build_section_geometry(
        station.layout,
        af,
        skin_thickness=0.005,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.006,
        outputs=("midlines",),
        web_x_positions_m=station.web_x_positions_m,
        web_y_coords=station.web_y_coords,
    )
    assert calls["n"] == 0


def test_shell_mesh_inputs_from_geometry_result_no_second_extract(monkeypatch):
    """Strips-first handoff must not re-run SDF midline export."""
    from blade_structure_utilities.mesh_utilities.section_mesh import handoff as mesh_mod
    from blade_structure_utilities.mesh_utilities.section_mesh.handoff import (
        _section_to_shell_inputs as from_sec_mod,
    )

    af = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    layout = resolve_system_type("1B-CN")
    result = build_section_geometry(
        layout,
        af,
        skin_thickness=0.005,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.006,
        twist_rad=0.1,
        outputs=("midlines",),
    )
    calls = {"n": 0}
    real = from_sec_mod.build_shell_midline_strips

    def _counting(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    monkeypatch.setattr(from_sec_mod, "build_shell_midline_strips", _counting)
    shell_inputs = mesh_mod.shell_mesh_inputs_from_geometry_result(result)
    assert calls["n"] == 0
    assert len(shell_inputs.midlines) == len(result.midline_strips)
    assert shell_inputs.chord_m == pytest.approx(float(result.chord_m))
    assert shell_inputs.twist_rad == pytest.approx(float(result.twist_rad))

    # Section→strips helper extracts once and matches strip count.
    from_section = from_sec_mod.shell_mesh_inputs_from_section(
        result.section, twist_rad=0.1
    )
    assert calls["n"] == 1
    assert len(from_section.midlines) == len(shell_inputs.midlines)


def test_cell_wall_loops_gate():
    af = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    layout = resolve_system_type("1B-CN")
    result = build_section_geometry(
        layout,
        af,
        skin_thickness=0.005,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.006,
        outputs=("cell_wall_loops",),
    )
    assert result.cell_wall_loops is not None
    assert len(result.cell_wall_loops) >= 1


def test_build_section_view_accepts_t_star_aliases():
    from blade_structure_utilities.geometry_utilities.section_geometry import (
        build_section_view,
    )

    af = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    layout = resolve_system_type("1B-CN")
    via_alias = build_section_view(
        af,
        layout,
        twist_angle_rad=0.0,
        skin_thickness=0.005,
        web_thickness=0.006,
        cap_height=0.008,
    )
    via_t = build_section_view(
        af,
        layout,
        twist_angle_rad=0.0,
        skin_thickness=0.005,
        shear_web_thickness=0.006,
        spar_cap_thickness=0.008,
    )
    assert "skin_suction_surface" in via_alias.labels
    assert via_alias.labels == via_t.labels
    with pytest.raises(ValueError, match="Conflicting thickness"):
        build_section_view(
            af,
            layout,
            twist_angle_rad=0.0,
            skin_thickness=0.005,
            shear_web_thickness=0.006,
            spar_cap_thickness=0.008,
            cap_height=0.009,
        )
