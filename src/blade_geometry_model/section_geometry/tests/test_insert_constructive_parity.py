"""LE/TE constructive emitters vs residual-phi occupancy (parity gate)."""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.constructive.emitters import (
    emit_leading_edge_insert_polygon,
    emit_trailing_edge_insert_polygon,
)
from blade_geometry_model.section_geometry.constructive.labels import (
    INSERTS_CONSTRUCTIVE_BY_DEFAULT,
    is_constructive_label,
    is_residual_label,
)
from blade_geometry_model.section_geometry.recipes.walls import (
    InsertWall,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)


def _shoelace_area(poly) -> float:
    if poly is None or getattr(poly, "is_empty", True):
        return 0.0
    return float(getattr(poly, "area", 0.0) or 0.0)


def _sdf_occupancy_area(comp, *, n: int = 220, pad: float = 0.02) -> float:
    af = getattr(comp, "_af", None)
    verts = np.asarray(af.vertices, dtype=float)
    xmin, ymin = verts.min(axis=0) - pad
    xmax, ymax = verts.max(axis=0) + pad
    xs = np.linspace(xmin, xmax, n)
    ys = np.linspace(ymin, ymax, n)
    xx, yy = np.meshgrid(xs, ys)
    phi = np.asarray(comp(xx.ravel(), yy.ravel()), dtype=float).reshape(xx.shape)
    dx = (xmax - xmin) / (n - 1)
    dy = (ymax - ymin) / (n - 1)
    return float(np.count_nonzero(phi <= 0.0)) * dx * dy


def test_inserts_constructive_by_default_flag():
    assert INSERTS_CONSTRUCTIVE_BY_DEFAULT is True
    assert is_constructive_label("leading_edge_insert")
    assert is_constructive_label("trailing_edge_insert")
    assert not is_residual_label("leading_edge_insert")


def test_insert_emitters_return_valid_or_empty_without_error():
    airfoil = AirfoilSDF.from_naca("2412", n_points=160, chord=1.0)
    te = InsertWall(
        kind="trailing_edge",
        skin_thickness=0.006,
        x_limit=0.85,
    )
    le = InsertWall(
        kind="leading_edge",
        skin_thickness=0.006,
        x_limit=0.08,
        radius=0.08,
    )
    te_poly = te.to_boundary_polygon(airfoil)
    le_poly = le.to_boundary_polygon(airfoil)
    assert te_poly is not None
    assert le_poly is not None
    assert _shoelace_area(te_poly) >= 0.0
    assert _shoelace_area(le_poly) >= 0.0


def test_insert_wall_allows_constructive_maturity():
    wall = InsertWall(
        kind="trailing_edge",
        skin_thickness=0.006,
        x_limit=0.8,
        emission_maturity="constructive",
    )
    assert wall.emission_maturity == "constructive"
    airfoil = AirfoilSDF.from_naca("0012", n_points=120, chord=1.0)
    poly = wall.to_boundary_polygon(airfoil)
    assert poly is not None


@pytest.mark.parametrize("layout_key", ["2D-F"])
def test_insert_constructive_vs_sdf_area_parity(layout_key):
    """Constructive insert area vs SDF occupancy within wall-solid hybrid band."""
    from blade_geometry_model.section_geometry.sections.multicell import (
        MultiCellSection,
    )

    airfoil = AirfoilSDF.from_naca("2412", n_points=180, chord=1.0)
    mcs = MultiCellSection.twin_web(
        airfoil,
        skin_thickness=0.006,
        web_thickness=0.008,
        cap_height=0.014,
        trailing_edge_insert_x_m=0.85,
        leading_edge_insert_x_m=0.08,
        le_radius=0.08,
        core_enabled=True,
    )
    comps = dict(getattr(mcs, "components_unrotated", None) or mcs.components)
    inserts = {
        k: v
        for k, v in comps.items()
        if "insert" in str(k).lower()
    }
    assert inserts, "expected LE/TE insert components"

    # Laminate hybrid area band (see DUAL_REP_ACCEPTANCE / hybrid gates).
    rtol, atol = 1e-3, 5e-4
    for key, comp in inserts.items():
        name = type(comp).__name__
        if name == "TEInsert":
            poly = emit_trailing_edge_insert_polygon(comp)
        elif name == "LEInsert":
            poly = emit_leading_edge_insert_polygon(comp)
        else:
            continue
        a_poly = _shoelace_area(poly)
        a_sdf = _sdf_occupancy_area(comp)
        if a_sdf < 1e-6 and a_poly < 1e-6:
            continue
        assert a_sdf >= 1e-6, f"{key}: SDF occupancy empty"
        abs_err = abs(a_poly - a_sdf)
        rel = abs_err / max(a_sdf, 1e-12)
        assert rel < rtol or abs_err < atol, (
            f"{key}: constructive area {a_poly:.6g} vs SDF {a_sdf:.6g} "
            f"(rel={rel:.3f}, abs={abs_err:.3g})"
        )
