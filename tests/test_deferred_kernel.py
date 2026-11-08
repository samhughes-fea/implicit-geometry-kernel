"""WallType leftovers, G9/WP6 polish, mid-surface CAD, and BladeConcept."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from implicit_geometry_kernel import (
    Airfoil,
    BladeConcept,
    BladeGeometry,
    SpanwiseScheme,
    Station,
    TransverseMember,
)
from implicit_geometry_kernel.extract import feature_curve
from implicit_geometry_kernel.sdf import LoftedPolygonSDF
from implicit_geometry_kernel.sdf.kernels.primitives import sdf_polygon


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations,
        "2D-F-LC-0R",
        0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48), pivot_fraction=0.3)


def test_named_feature_curve_matches_partition_edge():
    geo = _geo()
    g = geo.section_at(4.0).shell_graph()
    edges = g.partition_edges()
    assert edges
    e = edges[0]
    zs = np.linspace(0.0, 10.0, 5)
    named = feature_curve(geo, e.name_a, e.name_b, zs)
    from implicit_geometry_kernel.blade.shell_graph import partition_edge_curve

    direct = partition_edge_curve(geo, e.name_a, e.name_b, zs)
    np.testing.assert_allclose(named, direct)


def test_lofted_polygon_bvh_matches_unaccelerated_phi():
    geo = _geo()
    pts = geo.oml_polyline(4.0)
    query = np.column_stack([
        np.linspace(pts[:, 0].min() - 0.1, pts[:, 0].max() + 0.1, 40),
        np.linspace(pts[:, 1].min() - 0.05, pts[:, 1].max() + 0.05, 40),
        np.full(40, 4.0),
    ])
    fast = LoftedPolygonSDF(lambda z: geo.oml_polyline(float(z))[:, :2], accelerate=True)
    slow = LoftedPolygonSDF(lambda z: geo.oml_polyline(float(z))[:, :2], accelerate=False)
    np.testing.assert_allclose(fast.evaluate(query).phi, slow.evaluate(query).phi, atol=1e-12, rtol=1e-12)
    poly = geo.oml_polyline(4.0)[:, :2]
    np.testing.assert_allclose(
        slow.evaluate(query).phi,
        sdf_polygon(query[:, 0], query[:, 1], poly),
        atol=1e-12,
        rtol=1e-12,
    )


def test_write_mid_surfaces_obj_and_stl(tmp_path: Path):
    geo = _geo()
    obj = geo.write_mid_surfaces(tmp_path / "mids.obj", n_span=5, n_chord=12)
    stl = geo.write_mid_surfaces(tmp_path / "mids.stl", n_span=5, n_chord=12)
    text = Path(obj).read_text(encoding="utf-8")
    assert "g " in text and "v " in text and "f " in text
    assert "solid" in Path(stl).read_text(encoding="utf-8")


def test_write_mid_surfaces_step_skips_without_cadquery(tmp_path: Path):
    from implicit_geometry_kernel.blade.cad import cadquery_available

    geo = _geo()
    if not cadquery_available():
        with pytest.raises(ImportError, match="cadquery"):
            geo.write_mid_surfaces(tmp_path / "mids.step", n_span=3, n_chord=8)
        return
    path = geo.write_mid_surfaces(tmp_path / "mids.step", n_span=3, n_chord=8)
    assert Path(path).is_file()
    assert Path(path).stat().st_size > 0


def test_blade_concept_compiles_to_existing_kwargs():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    rib = TransverseMember("rib_mid", z=5.0, thickness=0.012)
    concept = BladeConcept("2D-F-LC-0R", SpanwiseScheme(transverse=(rib,)))
    blade = BladeGeometry.from_concept(
        concept,
        stations,
        0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48), pivot_fraction=0.3)
    legacy = BladeGeometry(
        stations,
        "2D-F-LC-0R",
        0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48),
        transverse=(rib,), pivot_fraction=0.3)
    assert blade.concept.layout == "2D-F-LC-0R"
    assert blade.concept.scheme.transverse == (rib,)
    assert len(blade.section_at(5.0).shell_graph().faces) == len(legacy.section_at(5.0).shell_graph().faces)
