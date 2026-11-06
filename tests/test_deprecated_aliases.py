"""Historical names stay until 0.6.0 and warn."""

from __future__ import annotations

import pytest


def test_cellloop_alias_warns() -> None:
    with pytest.warns(DeprecationWarning, match="CellLoop"):
        from implicit_geometry_kernel import CellLoop, CellPolygon

        assert CellLoop is CellPolygon


def test_spanweb_alias_warns() -> None:
    with pytest.warns(DeprecationWarning, match="SpanWeb"):
        from implicit_geometry_kernel import ShearWeb, SpanWeb

        assert SpanWeb is ShearWeb


def test_section_geometry_runner_alias_warns() -> None:
    with pytest.warns(DeprecationWarning, match="SectionGeometryRunner"):
        from implicit_geometry_kernel.section import SectionGeometry, SectionGeometryRunner

        assert SectionGeometryRunner is SectionGeometry


def test_station_at_alias_warns() -> None:
    from implicit_geometry_kernel.blade.airfoil import Airfoil
    from implicit_geometry_kernel.blade.loft import BladeGeometry
    from implicit_geometry_kernel.blade.stations import Station

    geo = BladeGeometry(
        [
            Station(z=0.0, chord=1.0, airfoil=Airfoil.naca4("0012"), twist=0.0),
            Station(z=1.0, chord=1.0, airfoil=Airfoil.naca4("0012"), twist=0.0),
        ],
        "0A",
        0.01,
    )
    with pytest.warns(DeprecationWarning, match="station_at"):
        a = geo.station_at(0.5)
    b = geo.interpolate(0.5)
    assert a[0] == b[0]


def test_web_x_positions_m_kwarg_warns() -> None:
    from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF
    from implicit_geometry_kernel.section.layout.system_layout import resolve_system_type
    from implicit_geometry_kernel.section.spec.section_spec import SectionSpec

    layout = resolve_system_type("2D-F")
    af = AirfoilSDF.from_naca("2412", n_points=40, chord=1.0)
    with pytest.warns(DeprecationWarning, match="web_x_positions_m"):
        spec = SectionSpec.from_kwargs(
            layout,
            af,
            skin_thickness=0.006,
            spar_cap_thickness=0.014,
            shear_web_thickness=0.008,
            web_x_positions_m=(0.25, 0.5),
        )
    assert spec.webs == (0.25, 0.5)
