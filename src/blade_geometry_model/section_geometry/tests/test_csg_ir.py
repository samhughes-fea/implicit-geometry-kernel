"""Regression tests for compiled CSG expression graphs."""

import numpy as np

from blade_geometry_model.section_geometry.sdf.ir import (
    Box,
    Circle,
    Intersect,
    Offset,
    Polygon,
    Shell,
    Slab,
    Subtract,
    Union,
    compile_expr,
    eval_expr,
    eval_expr_batch,
    intersect_expr,
    union_expr,
)
from blade_geometry_model.section_geometry.sdf.csg import intersect, offset, shell, subtract, union
from blade_geometry_model.section_geometry.sdf.grid import SDFGrid
from blade_geometry_model.section_geometry.sdf.primitives import sdf_box, sdf_circle
from blade_geometry_model.section_geometry import AirfoilSDF
from blade_geometry_model.section_geometry.recipes.sdf_ir import (
    LEInsert,
    OuterSkin,
    SandwichCore,
    ShearWeb,
    SparCap,
    TEInsert,
)


def test_compile_expr_matches_callable_csg():
    grid = SDFGrid.from_bbox(-2.0, 2.0, -2.0, 2.0, nx=120, ny=100)
    expr = Subtract(
        base=Union(Circle(cx=-0.2, cy=0.0, r=0.8), Box(cx=0.5, cy=0.0, half_w=0.5, half_h=0.4)),
        cutter=Offset(Circle(cx=0.0, cy=0.0, r=0.5), amount=-0.08),
    )
    phi_expr = compile_expr(expr, grid)

    f = subtract(
        union(
            lambda x, y: sdf_circle(x, y, cx=-0.2, cy=0.0, r=0.8),
            lambda x, y: sdf_box(x, y, cx=0.5, cy=0.0, half_w=0.5, half_h=0.4),
        ),
        offset(lambda x, y: sdf_circle(x, y, cx=0.0, cy=0.0, r=0.5), -0.08),
    )
    phi_ref = grid.eval(f)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_compile_expr_subtree_memoization_reuses_cache():
    grid = SDFGrid.from_bbox(-1.5, 1.5, -1.0, 1.0, nx=90, ny=70)
    leaf = Circle(cx=0.0, cy=0.0, r=0.6)
    expr = Intersect(Union(leaf, Box(cx=0.3, cy=0.0, half_w=0.4, half_h=0.3)), leaf)
    memo = {}
    phi = compile_expr(expr, grid, cache=memo)
    assert isinstance(phi, np.ndarray)
    assert leaf in memo
    assert expr in memo


def test_trailing_edge_insert_expr_matches_callable():
    af = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    comp = TEInsert(af, skin_thickness=0.003, x_start=0.72)
    assert comp._expr is not None
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=160, ny=90)
    phi_expr = compile_expr(comp._expr, grid)
    phi_ref = grid.eval(comp)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_leading_edge_insert_expr_matches_callable():
    af = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    comp = LEInsert(af, skin_thickness=0.003, x_end=0.12, le_x=0.0, le_y=0.0, radius=0.12)
    assert comp._expr is not None
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=160, ny=90)
    phi_expr = compile_expr(comp._expr, grid)
    phi_ref = grid.eval(comp)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_shear_web_expr_matches_callable_chord_normal():
    af = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    comp = ShearWeb(
        af,
        skin_thickness=0.003,
        x_top=0.25,
        y_top=0.06,
        x_bot=0.25,
        y_bot=-0.04,
        thickness=0.004,
        alignment="chord_normal",
        twist_angle=0.0,
    )
    assert comp._expr is not None
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=160, ny=90)
    phi_expr = compile_expr(comp._expr, grid)
    phi_ref = grid.eval(comp)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_shear_web_expr_matches_callable_flapwise():
    af = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    comp = ShearWeb(
        af,
        skin_thickness=0.003,
        x_top=0.25,
        y_top=0.06,
        x_bot=0.25,
        y_bot=-0.04,
        thickness=0.004,
        alignment="flapwise",
        twist_angle=0.2,
    )
    assert comp._expr is not None
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=160, ny=90)
    phi_expr = compile_expr(comp._expr, grid)
    phi_ref = grid.eval(comp)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_sandwich_core_expr_matches_callable_chord_clips():
    af = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    core = SandwichCore(
        af,
        skin_thickness=0.003,
        exclusion_sdfs=[],
        x_start=0.2,
        x_end=0.7,
    )
    assert core._expr is not None
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=160, ny=90)
    phi_expr = compile_expr(core._expr, grid)
    phi_ref = grid.eval(core)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_sandwich_core_expr_matches_callable_with_exclusion():
    af = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    web = ShearWeb(
        af,
        skin_thickness=0.003,
        x_top=0.3,
        y_top=0.05,
        x_bot=0.3,
        y_bot=-0.05,
        thickness=0.004,
    )
    core = SandwichCore(
        af,
        skin_thickness=0.003,
        exclusion_sdfs=[web],
        x_start=0.2,
        x_end=0.7,
    )
    assert core._expr is not None
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=160, ny=90)
    phi_expr = compile_expr(core._expr, grid)
    phi_ref = grid.eval(core)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_outer_skin_expr_matches_callable():
    af = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    comp = OuterSkin(af, thickness=0.006)
    assert comp.expr is not None
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=160, ny=90)
    phi_expr = compile_expr(comp.expr, grid)
    phi_ref = grid.eval(comp)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_spar_cap_expr_matches_callable():
    af = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    comp = SparCap(
        af,
        skin_thickness=0.003,
        x_start=0.2,
        x_end=0.6,
        cap_height=0.012,
        surface="upper",
    )
    assert comp._expr is not None
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=160, ny=90)
    phi_expr = compile_expr(comp._expr, grid)
    phi_ref = grid.eval(comp)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_continuous_spar_cap_expr_matches_callable_flapwise():
    af = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    comp = SparCap(
        af,
        skin_thickness=0.003,
        x_start=0.2,
        x_end=0.6,
        cap_height=0.012,
        surface="lower",
        twist_angle=0.15,
    )
    assert comp._expr is not None
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=160, ny=90)
    phi_expr = compile_expr(comp._expr, grid)
    phi_ref = grid.eval(comp)
    np.testing.assert_allclose(phi_expr, phi_ref, atol=1e-12, rtol=0.0)


def test_shell_and_slab_expr_match_callable():
    grid = SDFGrid.from_bbox(-1.5, 1.5, -1.0, 1.0, nx=90, ny=70)
    leaf = Circle(cx=0.0, cy=0.0, r=0.6)
    shell_expr = Shell(leaf, thickness=0.2)
    slab_expr = Slab(Offset(leaf, amount=-0.05), depth=0.08)
    phi_shell = compile_expr(shell_expr, grid)
    phi_slab = compile_expr(slab_expr, grid)
    f_shell = shell(lambda x, y: sdf_circle(x, y, cx=0.0, cy=0.0, r=0.6), 0.2)
    inner = offset(lambda x, y: sdf_circle(x, y, cx=0.0, cy=0.0, r=0.6), -0.05)

    def _slab(x, y):
        b = inner(x, y)
        return np.maximum(b, -b - 0.08)

    np.testing.assert_allclose(phi_shell, grid.eval(f_shell), atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(phi_slab, grid.eval(_slab), atol=1e-12, rtol=0.0)


def test_union_expr_intersect_expr_balanced():
    leaves = [Circle(cx=0.0, cy=0.0, r=0.3 + 0.05 * i) for i in range(5)]
    u = union_expr(leaves)
    i = intersect_expr(leaves)
    assert isinstance(u, Union)
    assert isinstance(i, Intersect)
    x = np.linspace(-1, 1, 40)
    y = np.linspace(-1, 1, 30)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    memo: dict = {}
    out_u = eval_expr_batch(u, xx, yy, memo)
    out_i = eval_expr_batch(i, xx, yy, memo)
    ref_u = leaves[0]
    for leaf in leaves[1:]:
        ref_u = Union(ref_u, leaf)
    ref_i = leaves[0]
    for leaf in leaves[1:]:
        ref_i = Intersect(ref_i, leaf)
    np.testing.assert_allclose(out_u, eval_expr(ref_u, xx, yy), atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(out_i, eval_expr(ref_i, xx, yy), atol=1e-12, rtol=0.0)


def test_eval_expr_batch_shares_polygon_subtree():
    verts = tuple((float(np.cos(t)), float(np.sin(t))) for t in np.linspace(0, 2 * np.pi, 32, endpoint=False))
    poly = Polygon(vertices=verts)
    inner = Offset(poly, amount=-0.02)
    expr = Intersect(Shell(poly, thickness=0.04), Slab(inner, depth=0.01))
    x = np.linspace(-1.2, 1.2, 50)
    y = np.linspace(-1.2, 1.2, 40)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    memo: dict = {}
    eval_expr_batch(expr, xx, yy, memo)
    assert poly in memo
    assert len([k for k in memo if isinstance(k, Polygon)]) == 1


def test_polygon_hash_is_precomputed():
    verts = tuple((float(i), float(i * i)) for i in range(200))
    poly = Polygon(vertices=verts)
    assert hash(poly) == poly._hash
    assert hash(poly) == hash(poly)

