"""Tests for shared section evaluation cache plumbing."""

import numpy as np

from blade_geometry_model.section_geometry.sdf.cache import SectionEvalCache
from blade_geometry_model.section_geometry.sdf.ir import Circle
from blade_geometry_model.section_geometry.sdf.grid import SDFGrid
from blade_geometry_model.section_geometry.props_emit.export import SectionPropertiesReport
from blade_geometry_model.section_geometry.sdf.transforms import rotate_field


class _CountingSDF:
    def __init__(self):
        self.calls = 0

    def __call__(self, x, y):
        self.calls += 1
        return np.sqrt((x - 0.5) ** 2 + y**2) - 0.2


class _DummySection:
    def __init__(self, comps):
        self._comps = comps

    def __iter__(self):
        return iter(self._comps)

    def __getitem__(self, label):
        return self._comps[label]


def test_section_properties_report_reuses_eval_cache():
    grid = SDFGrid.from_bbox(-1.0, 1.0, -1.0, 1.0, nx=64, ny=48)
    sdf = _CountingSDF()
    section = _DummySection({"component": sdf})
    cache = SectionEvalCache()

    SectionPropertiesReport(section, grid, backend="grid", eval_cache=cache)
    first = sdf.calls
    SectionPropertiesReport(section, grid, backend="grid", eval_cache=cache)
    second = sdf.calls

    assert first >= 1
    assert second == first


def test_eval_cache_owner_twist_fast_path_matches_rotated_callable():
    grid = SDFGrid.from_bbox(-1.0, 1.0, -1.0, 1.0, nx=80, ny=60)
    base = _CountingSDF()
    angle = 0.31
    rotated = rotate_field(base, angle)

    class _Owner:
        twist_rad = angle
        components_unrotated = {"component": base}

    cache = SectionEvalCache()
    phi_cached = cache.get_or_eval_with_owner("component", rotated, grid, owner=_Owner())
    phi_direct = grid.eval(rotated)
    np.testing.assert_allclose(phi_cached, phi_direct, atol=1e-12, rtol=0.0)


def test_eval_cache_uses_ir_by_default():
    grid = SDFGrid.from_bbox(-1.0, 1.0, -1.0, 1.0, nx=70, ny=50)

    class _Comp:
        def __init__(self):
            self.expr = Circle(cx=0.1, cy=-0.1, r=0.35)
            self._expr = self.expr

        def __call__(self, x, y):
            return np.sqrt((x - 0.1) ** 2 + (y + 0.1) ** 2) - 0.35

    class _Owner:
        twist_rad = 0.0
        components = {"component": _Comp()}

    cache = SectionEvalCache()
    phi_cached = cache.get_or_eval_with_owner(
        "component", _Owner.components["component"], grid, owner=_Owner()
    )
    phi_ref = grid.eval(_Owner.components["component"])
    np.testing.assert_allclose(phi_cached, phi_ref, atol=1e-12, rtol=0.0)

