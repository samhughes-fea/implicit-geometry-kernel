"""Wall recipe facade smoke tests."""

from __future__ import annotations

import pytest

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.constructive.emitters import (
    emit_sandwich_core_polygon,
)
from blade_geometry_model.section_geometry.constructive.export import (
    emit_constructive_polygons,
)
from blade_geometry_model.section_geometry.recipes import (
    SkinWall,
)
from blade_geometry_model.section_geometry.recipes.walls import (
    InsertWall,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)


def test_skin_wall_polygon_area_positive():
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    wall = SkinWall(thickness=0.006)
    poly = wall.to_boundary_polygon(airfoil)
    assert not poly.is_empty
    assert float(poly.area) > 0.0


def test_insert_wall_default_constructive_emits_boundary():
    airfoil = AirfoilSDF.from_naca("0012", n_points=100, chord=1.0)
    wall = InsertWall(kind="trailing_edge", skin_thickness=0.006, x_limit=0.75)
    assert wall.emission_maturity == "constructive"
    poly = wall.to_boundary_polygon(airfoil)
    assert poly is not None
    assert float(poly.area) >= 0.0


def test_insert_wall_residual_phi_defers_boundary():
    wall = InsertWall(
        kind="trailing_edge",
        skin_thickness=0.006,
        x_limit=0.75,
        emission_maturity="residual_phi",
    )
    with pytest.raises(NotImplementedError, match="constructive"):
        wall.to_boundary_polygon(None)


def test_insert_wall_constructive_emits_polygon():
    airfoil = AirfoilSDF.from_naca("0012", n_points=100, chord=1.0)
    wall = InsertWall(
        kind="trailing_edge",
        skin_thickness=0.006,
        x_limit=0.85,
        emission_maturity="constructive",
    )
    poly = wall.to_boundary_polygon(airfoil)
    assert poly is not None
    assert float(poly.area) >= 0.0


def test_skin_wall_to_sdf_expr_and_midline():
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    wall = SkinWall(thickness=0.006)
    expr = wall.to_sdf_expr(airfoil)
    assert expr is not None
    mid = wall.to_midline(airfoil)
    assert mid is not None and len(mid) >= 3


def test_core_wall_matches_emitter_area():
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    section = build_section_view(
        airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0
    )
    host = getattr(section, "host", section)
    recipes = dict(getattr(host, "wall_recipes", {}) or {})
    core_items = [
        (k, v) for k, v in recipes.items() if str(k).startswith("sandwich_core_")
    ]
    assert core_items
    polys = emit_constructive_polygons(section)
    comps = dict(getattr(host, "components_unrotated", None) or host.components)
    for key, wall in core_items:
        # Shared open cavity path uses exclusion already baked into export;
        # compare recipe↔direct emitter on the same SDF projection.
        comp = comps[key]
        via_recipe = wall.to_boundary_polygon(airfoil, wall_exclusion_union=None)
        via_emit = emit_sandwich_core_polygon(comp, None)
        assert abs(float(via_recipe.area) - float(via_emit.area)) <= 1e-12
        if key in polys and float(polys[key].area) > 1e-8:
            assert float(polys[key].area) > 0.0
