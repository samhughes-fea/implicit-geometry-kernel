"""Removed 0.6.0 names stay gone."""

from __future__ import annotations

import pytest

import implicit_geometry_kernel as k
from implicit_geometry_kernel import Airfoil
from implicit_geometry_kernel.blade.loft import BladeGeometry
from implicit_geometry_kernel.blade.stations import Station
from implicit_geometry_kernel.catalogue.system_layout import resolve_system_type
from implicit_geometry_kernel.section.spec.section_spec import SectionSpec


def test_cellloop_removed() -> None:
    with pytest.raises(ImportError):
        from implicit_geometry_kernel import CellLoop  # noqa: F401
    with pytest.raises(AttributeError):
        k.CellLoop


def test_spanweb_removed() -> None:
    with pytest.raises(ImportError):
        from implicit_geometry_kernel import SpanWeb  # noqa: F401
    with pytest.raises(AttributeError):
        k.SpanWeb


def test_section_geometry_runner_removed() -> None:
    import implicit_geometry_kernel.section as section

    with pytest.raises(ImportError):
        from implicit_geometry_kernel.section import SectionGeometryRunner  # noqa: F401
    with pytest.raises(AttributeError):
        section.SectionGeometryRunner


def test_station_at_removed() -> None:
    geo = BladeGeometry(
        [
            Station(z=0.0, chord=1.0, airfoil=Airfoil.naca4("0012"), twist=0.0),
            Station(z=1.0, chord=1.0, airfoil=Airfoil.naca4("0012"), twist=0.0),
        ],
        "0A-LC-0R",
        0.01, pivot_fraction=0.3)
    assert not hasattr(geo, "station_at")


def test_tag_node_wrappers_removed() -> None:
    import implicit_geometry_kernel.sdf.nodes.tagged as tagged

    for name in ("CircleSDF", "PolygonSDF", "ExprSDF"):
        assert not hasattr(tagged, name)


def test_wall_region_graph_cap_mode_removed() -> None:
    from implicit_geometry_kernel.blade.topology.graph import WallRegionGraph

    assert "cap_mode" not in WallRegionGraph.__dict__


def test_web_x_positions_m_kwarg_removed() -> None:
    layout = resolve_system_type("2D-F-LC-0R")
    af = Airfoil.naca4("2412").to_sdf(chord=1.0)
    with pytest.raises(TypeError, match="web_offsets"):
        SectionSpec.from_kwargs(
            layout,
            af,
            skin_thickness=0.006,
            spar_cap_thickness=0.014,
            shear_web_thickness=0.008,
            web_x_positions_m=(0.25, 0.5),
        )
