"""Ply-drop T(z,s), key points, OML export, queries, CSV/windIO, validate, nose, inserts."""

from pathlib import Path

import numpy as np
import pytest

from implicit_geometry_kernel import (
    Airfoil,
    BladeGeometry,
    PlyDrop,
    Station,
    read_planform_csv,
    write_windio,
)
from implicit_geometry_kernel.sdf.shapes.naca import naca_four_digit_thickness
from implicit_geometry_kernel.section.boundary.emitters import emit_trailing_edge_insert_polygon
from implicit_geometry_kernel.section.walls.walls import InsertWall
from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF


def _geo(**kw):
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    options = dict(spar_cap_thickness=0.008, shear_web_thickness=0.002, web_positions=(0.18, 0.48))
    options.update(kw)
    return BladeGeometry(stations, "2D-F", 0.003, **options)


def test_ply_drop_ramps_thickness_field():
    drop = PlyDrop("cap_drop", s0=0.2, s1=0.5, z0=2.0, z1=6.0, thickness=0.01, taper=0.2)
    geo = _geo(ply_drops=(drop,))
    assert float(geo.thickness_at(4.0, 0.35)) == pytest.approx(0.013)
    assert float(geo.thickness_at(4.0, 0.0)) == pytest.approx(0.003)
    assert float(geo.thickness_at(0.0, 0.35)) == pytest.approx(0.003)
    edge = float(geo.thickness_at(4.0, 0.2 - 0.1))
    assert 0.003 < edge < 0.013


def test_key_points_le_te_caps_webs():
    pts = _geo().key_points(5.0)
    assert "LE" in pts and "TE" in pts
    assert pts["TE"][0] > pts["LE"][0]
    assert any(k.startswith("cap_edge") for k in pts)
    assert any("shear_web" in k and "foot" in k for k in pts)


def test_write_oml_stl_and_obj(tmp_path: Path):
    geo = _geo()
    stl = geo.write_oml(tmp_path / "oml.stl", n_span=8, n_chord=16)
    obj = geo.write_oml(tmp_path / "oml.obj", n_span=8, n_chord=16)
    assert "facet" in Path(stl).read_text(encoding="utf-8")
    text = Path(obj).read_text(encoding="utf-8")
    assert text.startswith("v ") and "\nf " in text


def test_batched_query_and_validate():
    geo = _geo()
    cut = geo.section_at(5.0)
    p = cut.walls_of("skin")[0].points[3]
    out = geo.query([p, p + np.array([0.0, 0.0, 0.0])])
    assert out["thickness"].shape == (2,)
    assert out["distance"].shape == (2,)
    assert out["frame"].shape == (2, 3, 3)
    assert out["member"].shape == (2,)
    assert geo.validate() == []
    from implicit_geometry_kernel import TransverseMember

    broken = BladeGeometry(
        [Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
         Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3)],
        "2D-F", 0.003, transverse=(TransverseMember("rib", 50.0, 0.01),),
    )
    assert any("leaves the blade" in m for m in broken.validate())


def test_planform_csv_roundtrip_windio(tmp_path: Path):
    csv = tmp_path / "planform.csv"
    csv.write_text("z,chord,twist,pitch_axis,offset_x,offset_y,airfoil\n"
                   "0,2,0,0.3,0,0,0018\n"
                   "10,1,5,0.3,0,0.1,0012\n", encoding="utf-8")
    geo = read_planform_csv(csv, layout="2D-F", skin_thickness=0.004,
                            spar_cap_thickness=0.008, shear_web_thickness=0.004)
    assert geo.stations[0].chord == pytest.approx(2.0)
    doc = write_windio(geo)
    assert doc["outer_shape_bem"]["chord"]["values"][0] == pytest.approx(2.0)
    assert doc["outer_shape_bem"]["airfoil_position"]["labels"][0].startswith("NACA")


def test_naca_nose_thickness_matches_analytic():
    foil = Airfoil.naca4("0015", n=161)
    x = np.linspace(0.0, 0.01, 21)
    err = np.abs(foil.thickness(x) - naca_four_digit_thickness(x, 0.15))
    assert float(np.max(err)) < 1e-8


def test_insert_polygon_matches_shoelace():
    airfoil = AirfoilSDF.from_naca("2412", n_points=160, chord=1.0)
    wall = InsertWall(kind="trailing_edge", skin_thickness=0.006, x_limit=0.85)
    poly = wall.to_boundary_polygon(airfoil)
    assert poly is not None and not poly.is_empty
    xy = np.asarray(poly.exterior.coords, float)
    x, y = xy[:, 0], xy[:, 1]
    shoelace = 0.5 * float(np.abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))
    assert shoelace == pytest.approx(float(poly.area), abs=1e-12)
    emitted = emit_trailing_edge_insert_polygon(type("C", (), {
        "_af": airfoil, "skin_thickness": 0.006, "x_start": 0.85,
    })())
    assert float(emitted.area) == pytest.approx(float(poly.area), rel=1e-6)
