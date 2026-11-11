"""BladeConcept still compiles to BladeGeometry kwargs from blade.concepts."""

import numpy as np
import pytest

from implicit_geometry_kernel.blade.airfoil import Airfoil
from implicit_geometry_kernel.blade.concepts import BladeConcept, SpanwiseScheme
from implicit_geometry_kernel.blade.concepts import resolve_system_type
from implicit_geometry_kernel.blade.stations import Station, TransverseMember


def test_to_geometry_kwargs_includes_scheme():
    rib = TransverseMember("rib_mid", z=2.0, thickness=0.012)
    concept = BladeConcept("2D-F-LC-0R", SpanwiseScheme(transverse=(rib,)))
    kw = concept.to_geometry_kwargs()
    assert kw["layout"] == "2D-F-LC-0R"
    assert kw["transverse"] == (rib,)


def test_from_key_resolves_spec():
    concept = BladeConcept.from_key("2D-F-PA-10R")
    assert concept.layout.key == "2D-F-PA-10R"
    assert concept.layout.n_ribs == 10
    assert concept.layout.web_anchoring == "pitch_axis"
    assert resolve_system_type("2D-F-PA-10R").key == "2D-F-PA-10R"


def test_from_key_fills_ribs_when_stations_given():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 0.8, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    kw = BladeConcept.from_key("2D-F-LC-10R").to_geometry_kwargs(stations=stations)
    assert len(kw["transverse"]) == 10
    assert [m.name for m in kw["transverse"]] == [f"rib_{i}" for i in range(1, 11)]
    zs = [m.z for m in kw["transverse"]]
    assert zs == pytest.approx(np.linspace(0.0, 10.0, 12)[1:-1])
