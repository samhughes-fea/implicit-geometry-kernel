"""SS/PS skin solids and nested LE/bay/TE midline taxonomy."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
)
from blade_geometry_model.section_geometry.constructive import (
    emit_wall_solid_polygons,
    is_skin_label,
)
from blade_geometry_model.section_geometry.midlines.shell_midline_export import (
    build_shell_midline_strips,
)
from blade_geometry_model.section_geometry.junction_preparation.shell_midline_junction_preparation import (
    prepare_shell_midline_junctions,
)
from blade_geometry_model.section_geometry.cell_wall_loops.cell_wall_loops import (
    build_cell_wall_loops,
)
from blade_geometry_model.section_geometry.recipes.sdf_ir.outer_skin import (
    OuterSkin,
)


def test_is_skin_label_ss_ps():
    assert is_skin_label("outer_skin")
    assert is_skin_label("skin_suction_surface")
    assert is_skin_label("skin_pressure_surface")
    assert not is_skin_label("shear_web_0")


def test_ss_ps_skin_areas_sum_near_full_annulus():
    from blade_geometry_model.section_geometry.constructive.emitters import (
        emit_outer_skin_polygon,
    )

    af = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    t = 0.005
    full = OuterSkin(af, t)
    sec = MultiCellSection.twin_web(
        af,
        skin_thickness=t,
        web_thickness=0.006,
        cap_height=(0.014, 0.012),
        structural_family="D",
        core_enabled=False,
    )
    polys = emit_wall_solid_polygons(sec)
    a_ss = float(polys["skin_suction_surface"].area)
    a_ps = float(polys["skin_pressure_surface"].area)
    a_full = float(emit_outer_skin_polygon(full).area)
    assert a_ss > 0.0 and a_ps > 0.0
    # Camber clip + half-plane overlap near LE/TE can leave a small mismatch.
    np.testing.assert_allclose(a_ss + a_ps, a_full, rtol=0.08, atol=1e-4)


def test_twin_web_hierarchical_skin_arcs():
    af = AirfoilSDF.from_naca("2412", chord=1.0)
    sec = MultiCellSection.twin_web(
        af,
        skin_thickness=0.005,
        web_thickness=0.006,
        cap_height=(0.014, 0.012),
        structural_family="D",
        core_enabled=False,
    )
    raw = build_shell_midline_strips(sec, twist_rad=0.0)
    raw_skins = [s for s in raw if s.kind == "skin"]
    assert len(raw_skins) == 2
    # Asserted on the *raw* strips, which is now a guarantee the emit stage makes
    # rather than one `_ensure_ss_ps_skins` supplies a stage later. An SS/PS
    # half-ring is an open arc, so `closed` must be False here too â€” it used to be
    # True on both and was corrected downstream.
    assert {s.skin_side for s in raw_skins} == {"suction_surface", "pressure_surface"}
    assert all(not s.closed for s in raw_skins)

    topo = prepare_shell_midline_junctions(raw)
    skins = [s for s in topo if s.kind == "skin"]
    assert {s.skin_side for s in skins} == {"suction_surface", "pressure_surface"}
    # Twin-web: LE + bay_0 + TE per side â†’ 6 strips
    assert len(skins) == 6
    labels = {s.label for s in skins}
    assert "skin_suction_surface:leading_edge" in labels
    assert "skin_suction_surface:trailing_edge" in labels
    assert "skin_suction_surface:bay_0" in labels
    assert "skin_pressure_surface:leading_edge" in labels
    assert "skin_pressure_surface:trailing_edge" in labels
    assert "skin_pressure_surface:bay_0" in labels
    for s in skins:
        assert s.skin_side in ("suction_surface", "pressure_surface")
        assert s.skin_arc in ("leading_edge", "bay", "trailing_edge")
        if s.skin_arc == "bay":
            assert s.bay_index == 0
        else:
            assert s.bay_index is None

    loops = build_cell_wall_loops(topo)
    assert len(loops) == 3

