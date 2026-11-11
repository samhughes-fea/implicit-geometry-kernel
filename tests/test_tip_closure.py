"""Tip closure: every blade closes its skin at the tip, and it is not a rib."""

import warnings

import numpy as np
import pytest

from implicit_geometry_kernel import (
    Airfoil,
    BladeConcept,
    BladeGeometry,
    Station,
    TipClosure,
    TransverseMember,
    rib_face_polygons,
    tip_closure_polygons,
    topology_diff,
)
from implicit_geometry_kernel.blade.topology import (
    _template_from_cut,
    rib_cell_polygons_from_samples,
    sample_region_polylines,
)


def _stations():
    return [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]


def _geo(key="2D-F-LC-0R", skin=0.004, **kwargs):
    kw = dict(pivot_fraction=0.3, spar_cap_thickness=0.008, shear_web_thickness=0.004)
    kw.update(kwargs)
    return BladeGeometry(_stations(), key, skin, **kw)


@pytest.mark.parametrize("key", ["0A-LC-0R", "2D-F-LC-0R", "2D-F-LC-4R"])
def test_every_blade_has_a_tip_closure_in_the_skin_partition(key):
    geo = _geo(key)
    assert geo.tip_closure == TipClosure()
    assert geo.tip_closure_z == 10.0
    assert all(m.name != "tip_closure" for m in geo.transverse)
    faces = geo.shell_graph(10.0).faces
    tips = [f for f in faces if f.kind == "tip_closure"]
    assert tips
    for face in tips:
        assert face.surface == "skin"
        assert face.layers == ("skin",)
        assert face.thickness == pytest.approx(0.004)
        assert face.name.startswith("tip_closure")
    assert not [f for f in geo.shell_graph(9.0).faces if f.kind == "tip_closure"]


def test_tip_closure_takes_the_skin_thickness_law_at_the_tip():
    geo = _geo(skin=lambda z: 0.006 - 0.0002 * z)
    assert geo.tip_closure_thickness() == pytest.approx(0.004)
    face = next(f for f in geo.shell_graph(10.0).faces if f.kind == "tip_closure")
    assert face.thickness == pytest.approx(0.004)
    own = _geo(tip_closure=TipClosure(thickness=0.01))
    assert own.tip_closure_thickness() == pytest.approx(0.01)


def test_tip_closure_splits_at_webs_that_reach_the_tip():
    geo = _geo(web_positions=(0.18, 0.48))
    cut = geo.section_at(10.0)
    tips = [f for f in cut.shell_graph().faces if f.kind == "tip_closure"]
    assert len(tips) == len(cut.cells) == 3
    assert {f.name for f in tips} == {f"tip_closure:cell_{c.id}" for c in cut.cells}


def test_ribs_are_interior_and_the_tip_closure_is_extra():
    geo = _geo("2D-F-LC-10R")
    assert [m.name for m in geo.transverse] == [f"rib_{i}" for i in range(1, 11)]
    assert [m.z for m in geo.transverse] == pytest.approx(np.linspace(0.0, 10.0, 12)[1:-1])
    names = {n.split(":")[0] for n, _r in rib_face_polygons(geo, n_chord=8)}
    assert names == {f"rib_{i}" for i in range(1, 11)}
    assert {n.split(":")[0] for n, _r in tip_closure_polygons(geo, n_chord=8)} == {"tip_closure"}


def test_tip_closure_perimeter_is_the_skin_tip_ring():
    geo = _geo(web_positions=(0.18, 0.48))
    cut = geo.section_at(10.0)
    sg = cut.shell_graph()
    tmpl = _template_from_cut(geo.region_graph, cut, 10.0, 8)
    samples = sample_region_polylines(cut, sg, tmpl, geo.region_graph)
    rings = rib_cell_polygons_from_samples(tmpl, samples, "tip_closure")
    assert rings
    skin = {b.name for b in tmpl.loop}
    station = {tuple(np.round(p, 9)) for pts in samples.values() for p in np.asarray(pts, float)}
    ring_pts = {tuple(np.round(p, 9)) for _n, ring in rings for p in ring}
    # every closure vertex is a skin or web node of the tip station
    assert ring_pts <= station
    # and every skin node of the tip station is on the closure perimeter
    skin_pts = {tuple(np.round(p, 9)) for name, pts in samples.items() if name in skin for p in np.asarray(pts, float)}
    assert skin_pts <= ring_pts
    assert all(np.allclose(ring[:, 2], 10.0) for _n, ring in rings)


def test_tip_closure_is_a_topology_event_not_a_rib():
    geo = _geo()
    diff = topology_diff(geo.shell_graph(9.0), geo.shell_graph(10.0))
    kinds = {e.kind for e in diff.events if e.name.startswith("tip_closure")}
    assert kinds == {"gained_tip_closure"}


def test_concept_carries_the_tip_closure():
    concept = BladeConcept.from_key("2D-F-LC-0R")
    geo = BladeGeometry.from_concept(concept, _stations(), 0.004, pivot_fraction=0.3)
    assert geo.concept.scheme.tip_closure == geo.tip_closure
    tip = TipClosure(thickness=0.006)
    kw = BladeConcept("2D-F-LC-0R", type(concept.scheme)(tip_closure=tip)).to_geometry_kwargs()
    assert kw["tip_closure"] == tip


def test_legacy_tip_closure_member_becomes_the_tip_closure():
    with pytest.warns(DeprecationWarning, match="tip_closure"):
        geo = _geo(transverse=(TransverseMember("tip_closure", 10.0, 0.012),))
    assert geo.transverse == ()
    assert geo.tip_closure == TipClosure(thickness=0.012)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(ValueError, match="must sit at the tip"):
            _geo(transverse=(TransverseMember("tip_closure", 5.0, 0.012),))
        with pytest.raises(ValueError, match="not both"):
            _geo(transverse=(TransverseMember("tip_closure", 10.0, 0.012),), tip_closure=TipClosure())


def test_tip_closure_validates():
    with pytest.raises(ValueError):
        TipClosure(thickness=0.0)
    with pytest.raises(ValueError):
        TipClosure(name=" ")
