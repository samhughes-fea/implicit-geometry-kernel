"""Catalogue tokens: [W][X]-[Y]-[Z]-[R] with LC|LE|PA and required nR (including 0R)."""

import numpy as np
import pytest

from implicit_geometry_kernel import (
    SYSTEM_TYPE_KEYS,
    SYSTEM_TYPE_LAYOUT_KEYS,
    Airfoil,
    BladeGeometry,
    Station,
    resolve_system_type,
)
from implicit_geometry_kernel.catalogue.system_layout import (
    _LAYOUT_REGISTRY,
    effective_layout_for_web_count,
)


def _stations():
    return [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 0.8, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]


def test_bare_core_is_rejected():
    with pytest.raises(KeyError, match="Anchoring and rib count are required"):
        resolve_system_type("2D-F")


def test_missing_rib_count_is_rejected():
    with pytest.raises(KeyError, match="Rib count is required"):
        resolve_system_type("2D-F-LC")
    with pytest.raises(KeyError, match="Rib count is required"):
        resolve_system_type("2D-F-LC-R")


def test_every_token_round_trips():
    assert len(SYSTEM_TYPE_KEYS) == len(_LAYOUT_REGISTRY) * 3
    assert len(SYSTEM_TYPE_LAYOUT_KEYS) == len(_LAYOUT_REGISTRY)
    for key in SYSTEM_TYPE_KEYS:
        spec = resolve_system_type(key)
        assert spec.key == key
        assert key.endswith("-0R")
        assert spec.layout_key in _LAYOUT_REGISTRY
        assert spec.web_anchoring in ("local_chord", "leading_edge", "pitch_axis")
        assert spec.n_ribs == 0


def test_rib_count_round_trips():
    spec = resolve_system_type("2D-F-LC-10R")
    assert spec.key == "2D-F-LC-10R"
    assert spec.n_ribs == 10
    assert spec.web_anchoring == "local_chord"


def test_step_down_preserves_anchoring_and_ribs():
    spec = resolve_system_type("3D-F-PA-10R")
    stepped = effective_layout_for_web_count(spec, 2)
    assert stepped.key == "2D-F-PA-10R"
    assert stepped.web_anchoring == "pitch_axis"
    assert stepped.n_ribs == 10


def test_leading_edge_geometry_uses_offsets():
    stations = _stations()
    lc = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    le = BladeGeometry(stations, "2D-F-LE-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    assert le.web_offsets is not None
    assert lc.web_offsets is None
    assert le.web_offsets == pytest.approx(tuple(f * 2.0 for f in le.spec.web_chord_fracs))


def test_pitch_axis_holds_blade_x():
    stations = _stations()
    pa = BladeGeometry(stations, "2D-F-PA-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    assert pa._web_x_from_pitch is not None
    root = pa.section_at(0.0)
    tip = pa.section_at(5.0)
    def mean_x(cut, name):
        run = next(w for w in cut.shell_graph().webs if str(w.name) == name)
        return float(np.mean(np.asarray(run.points)[:, 0]))
    names = [w.name for w in root.shell_graph().webs]
    assert names
    for name in names:
        assert mean_x(tip, name) == pytest.approx(mean_x(root, name), abs=0.05 * 2.0)


def test_ribs_suffix_plants_transverse_members():
    stations = _stations()
    plain = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    ribbed = BladeGeometry(stations, "2D-F-LC-10R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    assert plain.transverse == ()
    assert len(ribbed.transverse) == 10
    assert all(m.name != "tip_closure" for m in ribbed.transverse)
    assert max(ribbed.forced_stations()) < 10.0
    assert plain.tip_closure == ribbed.tip_closure
