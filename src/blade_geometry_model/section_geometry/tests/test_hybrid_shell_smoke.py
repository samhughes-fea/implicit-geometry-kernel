"""Shell midline handoff still works when properties use the hybrid backend."""

from __future__ import annotations

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.midlines.shell_midline_export import (
    build_shell_midline_strips,
)
from blade_geometry_model.section_geometry.junction_preparation.shell_midline_junction_preparation import (
    ShellMidlineJunctionPreparationOptions,
    prepare_shell_midline_junctions,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF


def test_hybrid_props_then_shell_midlines():
    airfoil = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=6, max_cells=60000)
    props, meta = compute_section_properties(
        section, backend="hybrid", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    assert meta["backend"] == "hybrid"
    assert "skin_suction_surface" in props and "skin_pressure_surface" in props
    strips = build_shell_midline_strips(section, twist_rad=0.0)
    strips = prepare_shell_midline_junctions(strips, ShellMidlineJunctionPreparationOptions())
    assert len(strips) >= 3
    kinds = {s.kind for s in strips}
    assert "skin" in kinds
    skin_sides = {s.skin_side for s in strips if s.kind == "skin"}
    assert skin_sides == {"suction_surface", "pressure_surface"}
