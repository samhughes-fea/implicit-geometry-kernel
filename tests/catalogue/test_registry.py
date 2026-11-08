"""Catalogue tokens: [W][X]-[Y]-[Z]-[R] with LC|LE|PA and required nR."""

import pytest

from implicit_geometry_kernel.catalogue import (
    SYSTEM_TYPE_KEYS,
    SYSTEM_TYPE_LAYOUT_KEYS,
    effective_layout_for_web_count,
    resolve_system_type,
)
from implicit_geometry_kernel.catalogue.system_layout import _LAYOUT_REGISTRY


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


def test_section_layout_system_layout_shim_is_gone():
    import importlib
    from pathlib import Path

    import pytest

    pkg = Path(__file__).resolve().parents[2] / "src" / "implicit_geometry_kernel"
    assert not (pkg / "section" / "layout" / "system_layout.py").exists()
    with pytest.raises(ImportError):
        importlib.import_module("implicit_geometry_kernel.section.layout.system_layout")
