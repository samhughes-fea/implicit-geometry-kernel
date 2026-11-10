"""Characterization of parked Pass 2 gaps. Does not close them.

Numbers here are the inventory in ``docs/section/reports/MERGED_SHELL_GAPS.md``.
Loft from ``shell_graph``; snap with ``wall_mid_sdf`` (signed SDF of the
skin loop). Partition vertices lie on that loop.
"""

from pathlib import Path

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.topology.fields import StackThickness
from implicit_geometry_kernel.contract import ShellMidlineStrip
from implicit_geometry_kernel.extract import classify_trailing_edge
from implicit_geometry_kernel.sdf import LoftedPolygonSDF

ROOT = Path(__file__).resolve().parents[1]


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                         shear_web_thickness=0.002, web_positions=(0.18, 0.48), pivot_fraction=0.3)


def _catalogue():
    for path in (ROOT / "reference_blades", ROOT.parent / "reference-blades"):
        if (path / "wind" / "iea3.4_iea10" / "IEA-10-198-RWT.yaml").is_file() or (
            path / "iea3.4_iea10" / "IEA-10-198-RWT.yaml"
        ).is_file():
            return path
    return None


def _catalogue_yaml(rel: str):
    root = _catalogue()
    wind = root / "wind" / rel
    if wind.is_file():
        return wind
    return root / rel


def _max_abs_phi(field, pts):
    if len(pts) == 0:
        return 0.0
    return float(np.max(np.abs(field.evaluate(np.asarray(pts, float)).phi)))


def _loop_xy(loop):
    return np.vstack([r.points[:, :2] for r in loop])


# --- G1 ----------------------------------------------------------------------

def test_g1_solids_keep_caps_shell_adapter_does_not_duplicate():
    from implicit_geometry_kernel import shell_strips_from_graph

    cut = _geo().section_at(4.0)
    g = cut.shell_graph()
    strips = shell_strips_from_graph(cut)
    assert len(cut.walls) == 10
    assert len(g.runs()) == 8
    assert len(strips) == 8
    assert sum(1 for w in cut.walls if w.kind == "cap") == 2
    assert sum(1 for s in strips if s.kind == "cap") == 2
    assert {s.label for s in strips} == {r.name for r in g.runs()}
    assert "kind" in ShellMidlineStrip.__dataclass_fields__


# --- G2 ----------------------------------------------------------------------

def test_g2_from_meta_is_not_thickness_source():
    """D1.3: station_network_from_meta still carries thickness 0 (Wave 1)."""
    cut = _geo().section_at(4.0)
    from implicit_geometry_kernel.blade.topology.network import station_network_from_meta

    geo = _geo()
    field = geo.wall_meta()
    xs = np.concatenate([w.xy[:, 0] for w in cut.walls])
    ys = np.concatenate([w.xy[:, 1] for w in cut.walls])
    pad = 0.05
    bounds = (float(xs.min()) - pad, float(xs.max()) + pad, float(ys.min()) - pad, float(ys.max()) + pad)
    net = station_network_from_meta(geo.wall_mid_sdf(), field, 4.0, bounds, field.table, n=36)
    assert all(w.thickness == 0.0 for w in net.walls)
    assert cut.shell_graph().partition_edges()


# --- G3 ----------------------------------------------------------------------

def test_g3_loft_omits_missing_run_names():
    from implicit_geometry_kernel import shared_loft_names

    a = ("skin_bay", "shear_web_0", "cap_runout")
    b = ("skin_bay", "shear_web_0")
    shared = shared_loft_names(a, b)
    assert "cap_runout" not in shared
    assert shared == ("skin_bay", "shear_web_0")


# --- G4 ----------------------------------------------------------------------

def test_g4_stack_thickness_jumps_at_partition_edge():
    geo = _geo()
    cut = geo.section_at(4.0)
    g = cut.shell_graph()
    cap = next(r for r in g.loop if r.kind == "cap")
    skin = next(r for r in g.loop if r.kind == "skin")
    stack = StackThickness(geo)
    t_cap = float(stack(cap.points[len(cap.points) // 2].reshape(1, 3))[0])
    t_skin = float(stack(skin.points[len(skin.points) // 2].reshape(1, 3))[0])
    jump = abs(t_cap - t_skin)
    assert t_cap == pytest.approx(0.011)
    assert t_skin == pytest.approx(0.003)
    assert jump == pytest.approx(0.008)
    t_ends = stack(cap.points[[0, -1]])
    assert set(np.round(t_ends, 6)) == {0.003, 0.011}


# --- G5 ----------------------------------------------------------------------

def test_g5_sdf_vs_polyline_residual():
    geo = _geo()
    cut = geo.section_at(4.0)
    g = cut.shell_graph()
    mid = geo.wall_mid_sdf()
    cap_pts = np.vstack([r.points for r in g.loop if r.kind == "cap"])
    skin_pts = np.vstack([r.points for r in g.loop if r.kind == "skin"])
    cap_mid = next(r for r in g.loop if r.kind == "cap").points
    cap_mid_phi = abs(float(mid.evaluate(cap_mid[len(cap_mid) // 2].reshape(1, 3)).phi[0]))
    cap_phi = _max_abs_phi(mid, cap_pts)
    skin_phi = _max_abs_phi(mid, skin_pts)
    assert cap_mid_phi < 1e-12
    assert cap_phi < 1e-9
    assert skin_phi < 1e-9
    snapped = geo.snap_wall_mid(cap_pts)
    snap_phi = _max_abs_phi(mid, snapped.points)
    assert snap_phi < 1e-6
    assert snapped.converged.mean() > 0.9


# --- G6 ----------------------------------------------------------------------

def test_g6_web_foot_vs_stacked_loop():
    cut = _geo().section_at(4.0)
    g = cut.shell_graph()
    loop_xy = _loop_xy(g.loop)
    gaps = []
    skin_gaps = []
    for web in g.webs:
        for foot in (web.points[0, :2], web.points[-1, :2]):
            d_loop = float(np.min(np.linalg.norm(loop_xy - foot, axis=1)))
            gaps.append(d_loop)
            skins = np.vstack([w.xy for w in cut.walls if w.kind == "skin"])
            skin_gaps.append(float(np.min(np.linalg.norm(skins - foot, axis=1))))
    assert max(gaps) == pytest.approx(0.0)
    assert min(skin_gaps) == pytest.approx(0.0, abs=1e-9)
    assert len(gaps) == 4


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_g6_iea10_tip_section_at_builds_or_drops_web():
    import runpy

    ns = runpy.run_path(str(ROOT / "examples" / "blade" / "08_reference_blades" / "reference_blades.py"))
    _name, geo, _ = ns["load_windio"](_catalogue_yaml("iea3.4_iea10/IEA-10-198-RWT.yaml"))
    z_tip = float(geo.stations[-1].z)
    cut = geo.section_at(z_tip)
    assert cut.walls
    assert any(w.kind == "web" for w in cut.walls) or cut.dropped_webs


# --- G7 ----------------------------------------------------------------------

def test_g7_section_at_lru_and_lofted_aabb():
    from implicit_geometry_kernel.blade.topology.fields import CutWallMetaSDF

    geo = _geo()
    assert geo._section_cache_size == 256
    a = geo.section_at(4.0)
    b = geo.section_at(4.0 + 4e-10)
    assert a is b
    assert a.z == pytest.approx(4.0)
    meta = CutWallMetaSDF(geo)
    assert meta._cache_size == 256
    stack = StackThickness(geo)
    assert stack._cache_size == 256
    lofted = LoftedPolygonSDF(lambda z: geo.oml_polyline(z)[:, :2])
    assert lofted.accel_at(4.0) is not None
    assert lofted.bvh


# --- G8 ----------------------------------------------------------------------

def test_g8_example09_reads_windio_web_positions():
    src = (ROOT / "examples" / "blade" / "08_reference_blades" / "reference_blades.py").read_text(encoding="utf-8")
    assert "web_positions=(0.22, 0.50)" not in src
    assert "web_positions_fn_from_windio_structure" in src
    assert "cap_windows=(cap_windows if not web_positions else None)" not in src


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_g8_iea10_web_positions_from_yaml_anchors():
    import yaml

    from implicit_geometry_kernel.blade.windio import web_positions_from_windio_structure

    path = _catalogue_yaml("iea3.4_iea10/IEA-10-198-RWT.yaml")
    structure = yaml.safe_load(path.read_text(encoding="utf-8"))["components"]["blade"]["structure"]
    pos = web_positions_from_windio_structure(structure)
    assert pos is not None
    assert len(pos) == 3
    assert pos == tuple(sorted(pos))
    assert all(0.0 < f < 1.0 for f in pos)


# --- G9 ----------------------------------------------------------------------

def test_g9_cores_and_inserts_are_solids_not_loop_runs():
    cut = _geo().section_at(4.0)
    assert {w.kind for w in cut.walls} <= {"skin", "cap", "web"}
    g = cut.shell_graph()
    names = " ".join(r.name for r in g.runs())
    assert "core" not in names
    assert "insert" not in names
    assert {r.surface for r in g.runs()} <= {"skin", "web"}


def test_g9_opt_in_core_runs_from_polygon_boundaries():
    geo = _geo()
    g = geo.shell_graph(4.0, include_solids=True)
    cores = [r for r in g.solids if r.surface == "core"]
    assert cores
    assert all(len(r.points) >= 4 for r in cores)
    default = geo.shell_graph(4.0)
    assert default.solids == ()
    assert {r.surface for r in default.runs()} <= {"skin", "web"}


# --- G10 / G11 ----------------------------------------------------------------

def test_g10_thin_walled_not_in_this_repo():
    import importlib.util

    assert importlib.util.find_spec("implicit_geometry_kernel.thin_wall") is None


def test_g11_no_pcomp_in_kernel_src():
    src = ROOT / "src" / "implicit_geometry_kernel"
    hits = []
    for path in src.rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        if "PCOMP" in text or "pcomp" in text:
            hits.append(str(path.relative_to(ROOT)))
    assert hits == []


# --- WP6 leftovers ------------------------------------------------------------

def test_wp6_te_classifier_is_oml_not_merged_loop():
    geo = _geo()
    ring = geo.oml_polyline(4.0)
    oml = classify_trailing_edge(ring)
    g = geo.section_at(4.0).shell_graph()
    stacked = classify_trailing_edge(g)
    assert oml in {"sharp", "blunt"}
    assert stacked in {"sharp", "blunt"}
    assert oml == "sharp"
