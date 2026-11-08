"""Example viz helpers: suction-side names and twist-aware web fractions."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples"))

from _shared.blade_viz import is_suction_role, live_web_fracs, sample_planform  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.blade.topology import _chord_fracs  # noqa: E402


def test_is_suction_role_matches_authored_loop_names():
    assert is_suction_role("bay_su_0")
    assert is_suction_role("bay_su_1")
    assert is_suction_role("te_su")
    assert is_suction_role("le_su")
    assert is_suction_role("cap_su")
    assert is_suction_role("spar_cap_suction_surface")
    assert not is_suction_role("bay_pr_0")
    assert not is_suction_role("te_pr")
    assert not is_suction_role("fore_web")
    assert not is_suction_role("web_0")
    assert not is_suction_role("tip_closure:cell_0")
    assert not is_suction_role("rib_0")


def test_live_web_fracs_match_templates_along():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), twist=12.0, pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), twist=8.0, pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    cut = geo.section_at(5.0)
    live = live_web_fracs(cut)
    expected = []
    for run in cut.shell_graph().webs:
        pts = np.asarray(run.points, float)[:, :2]
        expected.append(float(np.clip(_chord_fracs(cut, pts[:1])[0], 0.0, 1.0)))
    assert live == expected
    assert live
    naive = []
    for run in cut.shell_graph().webs:
        xy = np.asarray(run.points, float)[:, :2]
        naive.append(float(np.mean(xy[:, 0] / max(cut.chord, 1e-9))))
    assert naive != live


def test_sample_planform_includes_rib_lines():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-4R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    data = sample_planform(geo, n=11)
    assert len(data["ribs"]) == 4
    zs = [row[0] for row in data["ribs"]]
    assert min(zs) > 0.0
    assert max(zs) == pytest.approx(10.0)
    for z, x_le, x_te, name in data["ribs"]:
        assert x_te > x_le
        assert name
