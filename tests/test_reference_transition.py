"""Published blades: cylinder–airfoil stations should cut after plateau collapse."""

import os
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def _catalogue():
    env = os.environ.get("REFERENCE_BLADES_DIR")
    roots = ([Path(env)] if env else []) + [ROOT / "reference_blades", ROOT.parent / "reference-blades"]
    for path in roots:
        if (path / "wind" / "nrel5mw" / "nrel5mw.yaml").is_file() or (
            path / "nrel5mw" / "nrel5mw.yaml"
        ).is_file():
            return path
    return None


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_nrel_transition_section_at():
    import runpy

    ns = runpy.run_path(str(ROOT / "examples" / "blade" / "08_reference_blades" / "reference_blades.py"))
    name, geo, _ = ns["load_windio"](
        ns["catalogue_yaml"](_catalogue(), "nrel5mw/nrel5mw.yaml")
    )
    cut = geo.section_at(3.07)
    assert len(cut.walls) >= 3
    kinds = {w.kind for w in cut.walls}
    assert "skin" in kinds


def _assert_unique_closed_loop(cut, needs_caps=True):
    g = cut.shell_graph()
    names = [r.name for r in g.loop]
    assert names, "empty skin loop"
    assert len(names) == len(set(names)), names
    assert g.closed_loop_chain() is not None, names
    scale = max(1.0, max(float(np.max(np.abs(r.points[:, :2]))) for r in g.loop))
    tol = 1e-6 * scale
    loop_xy = np.vstack([np.asarray(r.points[:, :2], float) for r in g.loop])
    for web in g.webs:
        for i in (0, -1):
            d = float(np.min(np.linalg.norm(loop_xy - web.points[i, :2], axis=1)))
            assert d <= tol, (web.name, i, d)
    if g.webs or needs_caps:
        assert any(r.kind == "cap" for r in g.loop)


_08_CASES = (
    "nrel5mw/nrel5mw.yaml",
    "iea15/IEA-15-240-RWT.yaml",
    "iea22/IEA-22-280-RWT.yaml",
    "iea3.4_iea10/IEA-3.4-130-RWT.yaml",
    "iea3.4_iea10/IEA-10-198-RWT.yaml",
)


def _08_ns():
    import runpy

    return runpy.run_path(str(ROOT / "examples" / "blade" / "08_reference_blades" / "reference_blades.py"))


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
@pytest.mark.parametrize(
    "rel, zs",
    [
        ("iea15/IEA-15-240-RWT.yaml", (0.0, 23.89)),
        ("iea3.4_iea10/IEA-3.4-130-RWT.yaml", (0.0, 16.28)),
    ],
)
def test_windio_loop_run_names_unique_at_root_and_max_chord(rel, zs):
    import runpy

    ns = runpy.run_path(str(ROOT / "examples" / "blade" / "08_reference_blades" / "reference_blades.py"))
    _name, geo, _ = ns["load_windio"](ns["catalogue_yaml"](_catalogue(), rel))
    for z in zs:
        _assert_unique_closed_loop(geo.section_at(float(z)))


def _max_chord_z(geo) -> float:
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = np.linspace(z0, z1, 81)
    chords = [geo.interpolate(float(z))[0] for z in zs]
    return float(zs[int(np.argmax(chords))])


def _painted_cap_span(cut) -> float:
    spans = []
    for run in cut.shell_graph().loop:
        if run.kind != "cap":
            continue
        xy = np.asarray(run.points[:, :2], float)
        spans.append(float(np.ptp(xy[:, 0])))
    assert spans, "no cap partition on the skin loop"
    return max(spans)


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
@pytest.mark.parametrize("z_key", ("root", "max_chord"))
def test_iea22_painted_cap_matches_yaml_window(z_key):
    import runpy

    from implicit_geometry_kernel.blade.windio import cap_window_table_from_windio_structure

    ns = runpy.run_path(str(ROOT / "examples" / "blade" / "08_reference_blades" / "reference_blades.py"))
    path = ns["catalogue_yaml"](_catalogue(), "iea22/IEA-22-280-RWT.yaml")
    _name, geo, z_of_s = ns["load_windio"](path)
    z = float(geo.stations[0].z) if z_key == "root" else _max_chord_z(geo)
    cut = geo.section_at(z)
    _assert_unique_closed_loop(cut)
    import yaml

    structure = yaml.safe_load(path.read_text(encoding="utf-8"))["components"]["blade"]["structure"]
    table = cap_window_table_from_windio_structure(structure, z_of_s)
    assert table is not None
    zt, vals = table
    lo = float(np.interp(z, zt, vals[:, 0]))
    hi = float(np.interp(z, zt, vals[:, 1]))
    expect = (hi - lo) * float(cut.chord)
    assert _painted_cap_span(cut) >= 0.5 * expect


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
@pytest.mark.parametrize("rel", _08_CASES)
def test_reference_blades_closed_conformal_loop_at_landmarks(rel):
    ns = _08_ns()
    _name, geo, _ = ns["load_windio"](ns["catalogue_yaml"](_catalogue(), rel))
    for _key, z, _label in ns["_landmarks"](geo):
        # caps bound by the webs exist only where webs are live; caps the file
        # places itself (cap windows) wherever it places them
        own = getattr(geo, "_cap_windows_fn", None)
        _assert_unique_closed_loop(geo.section_at(float(z)),
                                   needs_caps=own is not None and own(float(z)) is not None)


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_iea10_cap_loft_key_stable_along_span():
    ns = _08_ns()
    _name, geo, z_of_s = ns["load_windio"](
        ns["catalogue_yaml"](_catalogue(), "iea3.4_iea10/IEA-10-198-RWT.yaml")
    )
    zs = np.linspace(float(z_of_s[0, 1]), float(z_of_s[-1, 1]), 21)
    suction = set()
    pressure = set()
    for z in zs:
        try:
            g = geo.section_at(float(z)).shell_graph()
        except Exception:
            continue
        for r in g.loop:
            if r.kind != "cap":
                continue
            if "suction" in r.name:
                suction.add(r.name.split(":")[0])
            elif "pressure" in r.name:
                pressure.add(r.name.split(":")[0])
    assert suction == {"spar_cap_suction_surface"}
    assert pressure == {"spar_cap_pressure_surface"}


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_iea22_te_corners_do_not_chatter():
    from implicit_geometry_kernel.blade.topology import (
        ShellTopologyTemplate,
        _chord_fracs,
        _cut_tc,
        sample_region_polylines,
    )

    ns = _08_ns()
    _name, geo, _ = ns["load_windio"](ns["catalogue_yaml"](_catalogue(), "iea22/IEA-22-280-RWT.yaml"))
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    for z in np.linspace(z0, z1, 13):
        cut = geo.section_at(float(z))
        t = ShellTopologyTemplate.from_region_graph(
            geo.region_graph, float(z), n_chord=8, t_c=_cut_tc(cut),
        )
        samples = sample_region_polylines(cut, cut.shell_graph(), t, geo.region_graph)
        if "te_pr" not in samples:
            continue
        frac = float(_chord_fracs(cut, samples["te_pr"][-1:])[0])
        assert frac > 0.90, (z, frac)
        oml = np.asarray(geo.oml_polyline(float(z)), float)
        chord = max(float(cut.chord), 1e-6)
        hits = [np.asarray(samples["te_pr"][-1], float)]
        if "te" in samples:
            hits.append(np.asarray(samples["te"], float))
        stacked = np.vstack(hits)
        d0 = float(np.min(np.linalg.norm(stacked[:, :3] - oml[0, :3], axis=1)))
        d1 = float(np.min(np.linalg.norm(stacked[:, :3] - oml[-2, :3], axis=1)))
        assert min(d0, d1) < 0.08 * chord, (z, min(d0, d1), cut.chord)


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_iea10_cylinder_has_spar_webs_not_te_plug():
    from implicit_geometry_kernel.blade.topology import ShellTopologyTemplate

    ns = _08_ns()
    _name, geo, _ = ns["load_windio"](
        ns["catalogue_yaml"](_catalogue(), "iea3.4_iea10/IEA-10-198-RWT.yaml")
    )
    z0, span = float(geo.stations[0].z), float(geo.stations[-1].z) - float(geo.stations[0].z)
    names = []
    # the file's fore and rear webs start at 3 % span, the TE web at 5 %
    for z in (z0, z0 + 0.02 * span, z0 + 0.04 * span):
        cut = geo.section_at(float(z))
        from implicit_geometry_kernel.blade.topology import _cut_tc
        from implicit_geometry_kernel.blade.topology.ribs import _live_web_fracs

        t = ShellTopologyTemplate.from_region_graph(
            geo.region_graph, float(z), n_chord=8,
            live_web_fracs=_live_web_fracs(cut) or None, t_c=_cut_tc(cut),
        )
        names.append(tuple(w.name for w in t.webs))
    assert names[0] == () and names[1] == ()
    assert set(names[2]) == {"fore_web", "rear_web"}
    assert all("web_te" not in n or len(n) >= 2 for n in names)


def _span_metrics_for(geo, n=9):
    import sys

    sys.path.insert(0, str(ROOT / "examples"))
    from _shared.debug_mesh import as_quad, span_metrics
    from implicit_geometry_kernel.blade.topology import (
        loft_mapped_strips,
        sample_region_polylines,
        section_cut_near,
        templates_along,
    )

    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = np.linspace(z0, z1, n)
    tpls = templates_along(geo, zs, n_chord=8)
    walls = []
    for tmpl in tpls:
        cut = section_cut_near(geo, tmpl.z)
        if cut is None:
            continue
        samples = sample_region_polylines(cut, cut.shell_graph(), tmpl, geo.region_graph)
        walls.append((tmpl.z, tmpl, samples))
    faces, kinds = [], []
    for (_z0, t0, s0), (_z1, t1, s1) in zip(walls, walls[1:]):
        for name, face, _c in loft_mapped_strips(s0, s1, t0, t1):
            faces.append(face)
            kinds.append(name)
    face_arr = np.stack([as_quad(f) for f in faces]) if faces else np.zeros((0, 4, 3))
    return span_metrics(walls, face_arr, np.asarray(kinds, dtype="U64"), geo=geo)


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
@pytest.mark.parametrize(
    "rel, checks",
    [
        ("iea3.4_iea10/IEA-10-198-RWT.yaml", {
            # the file's webs start at 3 % span: none at the root station
            "n_web_clusters": 0, "n_coverage_holes": 0,
            "max_te_oml": 0.10, "max_te_chatter": 0.05,
        }),
        ("iea15/IEA-15-240-RWT.yaml", {"n_coverage_holes": 0}),
        ("iea22/IEA-22-280-RWT.yaml", {"n_coverage_holes": 0, "max_te_oml": 0.10}),
    ],
)
def test_windio_span_metrics(rel, checks):
    ns = _08_ns()
    _name, geo, _ = ns["load_windio"](ns["catalogue_yaml"](_catalogue(), rel))
    metrics = _span_metrics_for(geo)
    if "n_web_clusters" in checks:
        assert metrics["n_web_clusters"] == checks["n_web_clusters"], metrics
    if "max_wrap_frac" in checks:
        assert metrics["max_wrap_frac"] < checks["max_wrap_frac"], metrics
    if "n_coverage_holes" in checks:
        assert metrics["n_coverage_holes"] == checks["n_coverage_holes"], metrics
    if "max_te_oml" in checks:
        assert metrics["max_te_oml"] < checks["max_te_oml"], metrics
    if "max_te_chatter" in checks:
        assert metrics["max_te_chatter"] < checks["max_te_chatter"], metrics


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_iea10_loft_keepset_matches_intercept():
    from implicit_geometry_kernel.blade.topology import (
        sample_region_polylines,
        section_cut_near,
        templates_along,
    )

    ns = _08_ns()
    _name, geo, _ = ns["load_windio"](
        ns["catalogue_yaml"](_catalogue(), "iea3.4_iea10/IEA-10-198-RWT.yaml")
    )
    n_webs = len(geo.region_graph.webs)
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = np.unique(np.concatenate([np.linspace(z0, z1, 11), [20.0]]))
    tpls = templates_along(geo, zs, n_chord=8)
    root = tpls[0]
    assert not root.webs  # the file's webs start at 3 % span
    mid = next(t for t in tpls if abs(float(t.z) - 20.0) < 0.5)
    assert {w.index for w in mid.webs} == {0, 1, 2}
    assert {w.name for w in mid.webs} >= {"fore_web", "rear_web", "web_te"}
    for tmpl in tpls:
        cut = section_cut_near(geo, tmpl.z)
        dropped = {int(i) for i in cut.dropped_web_indices}
        expect = set(range(n_webs)) - dropped
        got = {int(w.index) for w in tmpl.webs}
        # a web is in the template where it is live (its span window) and not dropped
        assert got <= expect, (tmpl.z, got, expect, dropped)
        if 0.05 * (z1 - z0) <= float(tmpl.z) - z0 <= 0.475 * (z1 - z0):
            assert got == expect, (tmpl.z, got, expect, dropped)
        n_int = len(cut.shell_graph().webs)
        assert len(tmpl.webs) == n_int, (tmpl.z, len(tmpl.webs), n_int)
        samples = sample_region_polylines(cut, cut.shell_graph(), tmpl, geo.region_graph)
        for w in tmpl.webs:
            assert w.name in samples, (tmpl.z, w.name, sorted(samples))


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_iea34_fore_web_up_to_its_end():
    from implicit_geometry_kernel.blade.topology import (
        sample_region_polylines,
        section_cut_near,
        templates_along,
    )

    ns = _08_ns()
    _name, geo, _ = ns["load_windio"](
        ns["catalogue_yaml"](_catalogue(), "iea3.4_iea10/IEA-3.4-130-RWT.yaml")
    )
    # the file's webs end at 94 % span (z = 59.22 m)
    assert not templates_along(geo, (59.85,), n_chord=8)[0].webs
    zs = (58.5, 58.85, 59.2)
    tpls = templates_along(geo, zs, n_chord=8)
    assert len(tpls) == 3
    for tmpl in tpls:
        assert any(w.index == 0 or w.name == "fore_web" for w in tmpl.webs), tmpl.z
        cut = section_cut_near(geo, tmpl.z)
        samples = sample_region_polylines(cut, cut.shell_graph(), tmpl, geo.region_graph)
        names = {w.name for w in tmpl.webs}
        assert "fore_web" in samples or "web_0" in samples, (tmpl.z, names, sorted(samples))
