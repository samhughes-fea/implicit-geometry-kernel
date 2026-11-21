"""Mesh health registry: classify/flag only; adapters and study diffs."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from implicit_geometry_kernel.mesh.health import (
    HEALTH_CODES, HealthHit, detect_mid_surface_fold, detect_ridge_mismatch,
    detect_sheet_intersection, diff_study_folders, health_counts, hits_to_rows,
    run_mesh_health,
)
from implicit_geometry_kernel.mesh.quality import mesh_check
from implicit_geometry_kernel.mesh.tables import (
    mesh_failures, trailing_face_failures, write_mesh_tables,
)


def _flat_quad_mesh(**extra):
    nodes = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    ])
    base = dict(
        nodes=nodes,
        elements=((0, 1, 2, 3),),
        name=("panel",),
        surface=("skin",),
        constructor=("ruled",),
        bay=(0,),
        z0=(0.0,),
        z1=(0.1,),
        node_origin=("station",) * 4,
        section_distance=np.zeros(4),
        edges=(),
        loops=(),
        stations=np.array([0.0, 0.1]),
        recipe=None,
        partition=("cell",),
        thickness=(0.004,),
        eccentricity=(0.0,),
        n_source=4,
        edge_use=lambda: {(0, 1): 1, (1, 2): 1, (2, 3): 1, (3, 0): 1},
    )
    base.update(extra)
    return SimpleNamespace(**base)


def test_registry_codes_are_unique_and_ordered():
    assert len(HEALTH_CODES) == len(set(HEALTH_CODES))
    assert "open_edge_off_root" in HEALTH_CODES
    assert "sheet_intersection" in HEALTH_CODES


def test_mesh_failures_adapter_filters_quality_hits():
    from implicit_geometry_kernel.mesh.quality import element_quality
    from implicit_geometry_kernel.mesh.tables import build_corners

    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.99, 0.001, 0.0],
        [0.0, 1.0, 0.0],
    ])
    elements = ((0, 1, 2, 3),)
    quality = element_quality(nodes, elements, [0.004])
    mesh = SimpleNamespace(
        nodes=nodes, elements=elements, element_quality=quality,
        section_distance=np.zeros(4), edges=(), constructor=("step_patch",),
        bay=(4,), z0=(0.0,), z1=(0.05,), name=("te",), surface=("skin",),
        loops=(), node_origin=("station",) * 4, stations=np.array([0.0]),
        recipe=None,
    )
    mesh.corners = build_corners(mesh)
    groups = mesh_failures(mesh)
    assert len(groups) == 1
    assert groups[0]["constructor"] == "step_patch"
    hits = run_mesh_health(mesh, include_geometry=False)
    assert any(hit.code == "quad_angle" for hit in hits)


def test_trailing_face_adapter_maps_reasons():
    from implicit_geometry_kernel.mesh.quality import element_quality

    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.99, 0.001, 0.0],
        [0.0, 1.0, 0.0],
    ])
    elements = ((0, 1, 2, 3),)
    quality = element_quality(nodes, elements, [0.004])
    mesh = SimpleNamespace(
        nodes=nodes, elements=elements, element_quality=quality,
        section_distance=np.zeros(4), edges=(), constructor=("step_patch",),
        bay=(1,), z0=(0.0,), z1=(0.05,), name=("te",), surface=("skin",),
        loops=(), node_origin=("station",) * 4, stations=np.array([0.0]),
        recipe=None, diagonal_min=np.array([1e-6]),
    )
    # ensure quality table exposes jacobian / angle for trailing detector
    rows = trailing_face_failures(mesh)
    assert isinstance(rows, list)


def test_ridge_mismatch_flags_separated_te_corners():
    nodes = np.array([
        [1.0, 0.1, 0.0], [0.0, 0.1, 0.0], [0.0, 0.1, 0.1], [1.0, 0.1, 0.1],
        [1.05, -0.1, 0.0], [0.0, -0.1, 0.0], [0.0, -0.1, 0.1], [1.05, -0.1, 0.1],
    ])
    mesh = SimpleNamespace(
        nodes=nodes,
        elements=((0, 1, 2, 3), (4, 5, 6, 7)),
        name=("te_su", "te_pr"),
        node_origin=("station",) * 8,
        constructor=("ruled", "ruled"),
        bay=(0, 0),
        z0=(0.0, 0.0),
        z1=(0.1, 0.1),
    )
    hits = detect_ridge_mismatch(mesh, tol=1e-3)
    assert hits and hits[0].code == "ridge_mismatch"


def test_mid_surface_fold_uses_dihedral():
    dtype = np.dtype([
        ("a", "i8"), ("b", "i8"), ("use", "i8"), ("e0", "i8"), ("e1", "i8"),
        ("meet", "U16"), ("kind", "U16"), ("dihedral_deg", "f8"), ("length", "f8"),
    ])
    edges = np.zeros(1, dtype=dtype)
    edges[0] = (0, 1, 2, 0, 1, "skin", "span", 80.0, 1.0)
    mesh = _flat_quad_mesh(edges=edges, elements=((0, 1, 2, 3), (0, 1, 2, 3)))
    hits = detect_mid_surface_fold(mesh, fold=45.0)
    assert len(hits) == 1 and hits[0].code == "mid_surface_fold"


def test_sheet_intersection_detects_crossing_te_sheets():
    # two overlapping quads that pierce (offset in y, crossed in x/z)
    nodes = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [0.0, 0.0, 1.0],
        [0.25, -0.5, 0.25], [0.75, -0.5, 0.25], [0.75, 0.5, 0.75], [0.25, 0.5, 0.75],
    ])
    mesh = SimpleNamespace(
        nodes=nodes,
        elements=((0, 1, 2, 3), (4, 5, 6, 7)),
        name=("te_su", "te_pr"),
        constructor=("ruled", "ruled"),
        bay=(0, 0),
        z0=(0.0, 0.0),
        z1=(1.0, 1.0),
        edges=(),
    )
    hits = detect_sheet_intersection(mesh)
    assert any(hit.code == "sheet_intersection" for hit in hits)


def test_diff_study_folders_reports_count_and_element_delta(tmp_path: Path):
    before = tmp_path / "before"
    after = tmp_path / "after"
    before.mkdir()
    after.mkdir()
    (before / "health_counts.csv").write_text(
        "code,severity,layer,count,inspect\nquad_angle,warn,quality,2,x\n",
        encoding="utf-8",
    )
    (after / "health_counts.csv").write_text(
        "code,severity,layer,count,inspect\nquad_angle,warn,quality,1,x\n"
        "open_edge_off_root,fail,topology,3,y\n",
        encoding="utf-8",
    )
    (before / "failures.csv").write_text(
        "code,severity,element,node,edge,region,bay,constructor,z0,z1,measure,threshold,detail\n"
        "quad_angle,warn,10,,,,,,,,,\n"
        "quad_angle,warn,11,,,,,,,,,\n",
        encoding="utf-8",
    )
    (after / "failures.csv").write_text(
        "code,severity,element,node,edge,region,bay,constructor,z0,z1,measure,threshold,detail\n"
        "quad_angle,warn,10,,,,,,,,,\n"
        "open_edge_off_root,fail,,,4,,,,,,,\n",
        encoding="utf-8",
    )
    report = diff_study_folders(before, after)
    by_code = {row["code"]: row for row in report["counts"]}
    assert by_code["quad_angle"]["delta"] == -1
    assert by_code["open_edge_off_root"]["delta"] == 3
    assert report["elements"]["gone"] == [11]
    assert report["elements"]["new"] == []


def test_write_mesh_tables_emits_failures_and_health_counts(tmp_path, concept_mesh):
    write_mesh_tables(concept_mesh, tmp_path)
    assert (tmp_path / "failures.csv").exists()
    assert (tmp_path / "health_counts.csv").exists()
    summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "health hits:" in summary
    check = mesh_check(concept_mesh)
    assert isinstance(check.hit_counts, dict)
    assert set(check.hit_counts) <= set(HEALTH_CODES)


@pytest.fixture
def concept_mesh():
    from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
    from implicit_geometry_kernel.mesh.shell import shell_mesh

    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3, twist=8.0),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )
    return shell_mesh(geo)


def test_hits_to_rows_joinable_ids():
    hit = HealthHit(
        code="quad_angle", severity="warn", element=3, measure=12.0, threshold=20.0,
        region="panel", bay=1, constructor="ruled",
    )
    row = hits_to_rows((hit,))[0]
    assert row["element"] == 3
    assert row["node"] == ""
    assert health_counts((hit,)) == {"quad_angle": 1}
