"""Quadtree area mesh: circle/airfoil samplers, refinement, contour closure."""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.area_2d import (
    QuadtreeSpec,
    build_area_mesh_2d,
    build_area_meshes_2d,
    contour_from_quadtree_marching_squares,
    refine_quadtree,
    build_root_quadtree,
)


def _circle_sampler(radius: float = 0.2):
    def sampler(pts: np.ndarray) -> np.ndarray:
        arr = np.asarray(pts, dtype=float).reshape(-1, 2)
        return np.hypot(arr[:, 0], arr[:, 1]) - radius

    return sampler


def test_quadtree_refines_near_zero_contour():
    spec = QuadtreeSpec(-0.5, 0.5, -0.5, 0.5, max_depth=5, min_cell_size=0.02)
    mesh = build_root_quadtree(spec, _circle_sampler())
    mesh = refine_quadtree(mesh, _circle_sampler())
    leaves = mesh.leaf_cells()
    assert len(leaves) > 1
    min_w = min(c.width for c in leaves)
    assert min_w < 0.25


def test_contour_from_circle_is_closed():
    spec = QuadtreeSpec(-0.5, 0.5, -0.5, 0.5, max_depth=6, min_cell_size=0.01)
    field = build_area_mesh_2d(_circle_sampler(), spec)
    assert field.contour is not None
    c = field.contour
    assert np.allclose(c[0], c[-1], atol=1e-6)
    radii = np.hypot(c[:-1, 0], c[:-1, 1])
    assert np.allclose(radii, 0.2, atol=0.02)


def test_quadtree_spec_from_points():
    pts = np.array([[0.0, -0.05], [1.0, 0.05], [0.5, 0.12]])
    spec = QuadtreeSpec.from_points(pts, padding=0.1, min_cell_size=0.01)
    assert spec.x_min < 0.0
    assert spec.x_max > 1.0


def test_build_area_meshes_2d_multiple_fields():
    fields = {"a": _circle_sampler(0.15), "b": _circle_sampler(0.25)}
    spec = QuadtreeSpec(-0.5, 0.5, -0.5, 0.5, max_depth=4, min_cell_size=0.05)
    area = build_area_meshes_2d(fields, spec)
    assert set(area.fields.keys()) == {"a", "b"}
    assert area.fields["a"].stats["leaf_count"] >= 1


@pytest.mark.parametrize("code", ["0012"])
def test_airfoil_component_quadtree_contour(code: str):
    from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
    from blade_geometry_model.section_geometry.api.factory import (
        build_section_view,
    )
    from blade_geometry_model.section_geometry.layout.system_layout import (
        resolve_system_type,
    )
    from blade_geometry_model.section_geometry.props_emit.area_mesh import (
        build_section_area_meshes,
    )

    airfoil = AirfoilSDF.from_naca(code, n_points=80, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-F"), twist_angle_rad=0.0)
    verts = np.asarray(airfoil.vertices, dtype=float)
    spec = QuadtreeSpec.from_points(verts, padding=0.08, min_cell_size=0.012, max_depth=6)

    meshes = build_section_area_meshes(section, spec, twist_rad=0.0)
    skin = meshes.fields.get("skin_outer") or meshes.fields.get("outer_skin")
    if skin is None:
        skin = next(iter(meshes.fields.values()))
    if skin.contour is not None:
        assert np.allclose(skin.contour[0], skin.contour[-1], atol=1e-5)
