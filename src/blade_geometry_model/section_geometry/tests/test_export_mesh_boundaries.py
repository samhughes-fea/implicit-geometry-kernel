"""Export component boundaries from shared quadtree mesh (Phase 9)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.props_emit.export import (
    SectionPropertiesReport,
    export_section_json,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.sdf.grid import SDFGrid


def test_export_reuses_properties_report_and_mesh_boundaries(tmp_path: Path):
    airfoil = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=5, max_cells=50000)
    props, meta, mesh = compute_section_properties(
        section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=airfoil, with_mesh=True
    )
    assert mesh is not None
    report = SectionPropertiesReport.from_precomputed(section, props, meta, mesh)
    grid = SDFGrid.from_bbox(spec.x_min, spec.x_max, spec.y_min, spec.y_max, 48, 24)
    out = tmp_path / "section.json"
    export_section_json(
        section,
        grid,
        {},
        str(out),
        include_component_zero_contours=True,
        include_geometry_detail=True,
        properties_report=report,
    )
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert "section_properties" in payload
    assert "geometry" in payload
    n_with_boundary = 0
    for label in section.labels:
        lab = str(getattr(label, "value", label))
        boundary = payload["geometry"]["components"][lab]["boundary"]
        assert isinstance(boundary, list)
        if boundary:
            n_with_boundary += 1
            first = boundary[0]
            assert isinstance(first, list) and len(first) >= 2
            assert len(first[0]) == 2
    # At least the primary solids should export a boundary (web_0 can be empty at coarse depth)
    assert n_with_boundary >= 3
    assert payload["geometry"]["components"]["skin_suction_surface"]["boundary"]
    assert payload["geometry"]["components"]["skin_pressure_surface"]["boundary"]
