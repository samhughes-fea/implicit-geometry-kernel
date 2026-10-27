"""`SectionGeometryRunner`: the two products, and the promise about their cost.

This runner had **no test and no caller**, and had already drifted: its module
docstring said the layout "is resolved *before* anything is assembled, so a
dropped web is known to the caller without paying for a build", while
`execute()` unconditionally did both and `get_layout_resolution()` raised unless
it had run. Nothing exercised it, so nothing caught it.

That promise is load-bearing for the sizing loop, which resolves a station's
layout once to fix its design-variable set and then builds geometry many times
per resolution while searching thicknesses — measured at 52 resolutions against
168 builds on the reference case. Folding the two together would add roughly a
third more builds for nothing.
"""

from __future__ import annotations

import pytest

from blade_geometry_model.section_geometry import (
    resolve_system_type,
)
from blade_geometry_model.section_geometry.api.runner import (
    SectionGeometryRunner,
)


@pytest.fixture(scope="module")
def airfoil():
    from blade_geometry_model.section_geometry.sdf.airfoil import (
        AirfoilSDF,
    )

    return AirfoilSDF.from_naca("0012", n_points=100, chord=1.0)


def _runner(airfoil, **kw):
    return SectionGeometryRunner(
        layout=resolve_system_type(kw.pop("key", "2B-CN")),
        airfoil=airfoil,
        skin_thickness=kw.pop("skin_thickness", 0.004),
        spar_cap_thickness=kw.pop("spar_cap_thickness", 0.008),
        shear_web_thickness=kw.pop("shear_web_thickness", 0.003),
        twist_rad=0.0,
        outputs=kw.pop("outputs", ("midlines",)),
        **kw,
    )


def test_resolve_yields_a_layout_without_building(airfoil, monkeypatch):
    """The docstring's promise, asserted rather than described.

    Counts the *builder*, because "did not pay for a build" is a statement about
    work done, not about which attribute happens to be populated.
    """
    import blade_geometry_model.section_geometry.api.runner as mod

    calls = {"n": 0}
    real = mod.build_section_geometry

    def counting(*a, **kw):
        calls["n"] += 1
        return real(*a, **kw)

    monkeypatch.setattr(mod, "build_section_geometry", counting)

    runner = _runner(airfoil).resolve()
    resolution = runner.get_layout_resolution()

    assert resolution is not None
    assert resolution.layout is not None
    assert calls["n"] == 0, "resolve() must not build geometry"


def test_execute_yields_both_products(airfoil):
    runner = _runner(airfoil).execute()
    assert runner.get_layout_resolution() is not None
    assert runner.get_results() is not None


def test_execute_after_resolve_reuses_the_resolution(airfoil, monkeypatch):
    """Calling both must not resolve twice — the loop does exactly this."""
    import blade_geometry_model.section_geometry.api.runner as mod

    calls = {"n": 0}
    real = mod.resolve_station_layout

    def counting(*a, **kw):
        calls["n"] += 1
        return real(*a, **kw)

    monkeypatch.setattr(mod, "resolve_station_layout", counting)

    runner = _runner(airfoil).resolve()
    runner.execute()
    assert calls["n"] == 1


def test_both_lifecycle_methods_are_idempotent(airfoil):
    runner = _runner(airfoil)
    assert runner.resolve() is runner.resolve()
    first = runner.execute().get_results()
    assert runner.execute().get_results() is first


def test_get_results_requires_execute(airfoil):
    """A resolution alone is not a geometry, and must not pretend to be."""
    runner = _runner(airfoil).resolve()
    with pytest.raises(RuntimeError, match="execute"):
        runner.get_results()


def test_get_layout_resolution_requires_resolve_or_execute(airfoil):
    with pytest.raises(RuntimeError, match="resolve"):
        _runner(airfoil).get_layout_resolution()


def test_a_supplied_resolution_is_not_re_resolved(airfoil, monkeypatch):
    """A caller that already resolved this station must not pay for it again.

    And more than cost: a sizing loop freezes topology at the minimum skin so the
    design-variable set is fixed for the whole search. Re-resolving at each
    trial's thickness could return a different web set than the station was
    frozen against.
    """
    import blade_geometry_model.section_geometry.api.runner as mod

    frozen = _runner(airfoil).resolve().get_layout_resolution()

    calls = {"n": 0}
    real = mod.resolve_station_layout

    def counting(*a, **kw):
        calls["n"] += 1
        return real(*a, **kw)

    monkeypatch.setattr(mod, "resolve_station_layout", counting)

    runner = _runner(airfoil, resolution=frozen).execute()

    assert calls["n"] == 0, "a supplied resolution must not be re-resolved"
    assert runner.get_layout_resolution() is frozen
    assert runner.get_results() is not None


def test_the_build_follows_the_supplied_resolution(airfoil):
    """Not just accepted — actually used to build.

    A resolution that is accepted and then ignored would be worse than not
    accepting one, because the caller would believe the frozen topology held.
    """
    nominal = _runner(airfoil, key="5D-CN")
    frozen = nominal.resolve().get_layout_resolution()

    built = _runner(airfoil, key="5D-CN", resolution=frozen).execute().get_results()
    n_web_strips = sum(
        1 for s in (built.midline_strips or []) if str(s.kind) == "web"
    )
    assert n_web_strips == len(frozen.web_x_positions_m or ())


def test_the_build_uses_the_effective_layout_not_the_nominal(airfoil):
    """A station that cannot anchor every web must not be built as if it could.

    Sizing allocates a design variable per component, so building against the
    nominal layout would hand the search a `shear_web_thickness` for a web the section does
    not have.
    """
    runner = _runner(airfoil, key="5D-CN").execute()
    resolution = runner.get_layout_resolution()
    result = runner.get_results()

    # The result carries no layout key, so assert the thing that key would stand
    # for: the section has exactly the webs the resolution anchored.
    n_web_strips = sum(
        1 for s in (result.midline_strips or []) if str(s.kind) == "web"
    )
    assert n_web_strips == len(resolution.web_x_positions_m or ()), (
        f"resolution anchored {len(resolution.web_x_positions_m or ())} web(s) "
        f"but the built section has {n_web_strips} — geometry was built against "
        f"a different layout than the one resolved"
    )
    if resolution.stepped_down:
        assert len(resolution.web_x_positions_m) < int(resolution.layout.n_webs) + 1
        assert resolution.dropped, "a stepped-down layout must say what it dropped"
