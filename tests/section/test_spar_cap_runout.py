"""Spar-cap run-out gate: a cap that cannot be built at gauge is terminated."""

from __future__ import annotations

import math
from dataclasses import dataclass

import pytest

from implicit_geometry_kernel.catalogue.system_layout import resolve_system_type
from implicit_geometry_kernel.sdf import AirfoilSDF
from implicit_geometry_kernel.section.layout.spar_cap_runout import (
    apply_spar_cap_runout,
    max_spar_from_cavity,
)


@dataclass(frozen=True)
class _Resolution:
    layout: object
    web_x_positions_m: tuple[float, ...] = ()


def _naca(code: str, chord: float = 1.0) -> AirfoilSDF:
    return AirfoilSDF.from_naca(code, n_points=200, chord=chord)


def test_cavity_bound_keeps_the_clamp_void() -> None:
    """SG-039: d = h_local - t_s; the void is max(gap_min, 0.15 d)."""
    # d = 96 mm: the 15 % void (14.4 mm) governs.
    assert max_spar_from_cavity(0.10, 0.004, gap_min=0.002) == pytest.approx(0.5 * (0.096 - 0.0144))
    # d = 6 mm: the floor governs.
    assert max_spar_from_cavity(0.010, 0.004, gap_min=0.002) == pytest.approx(0.5 * (0.006 - 0.002))
    assert max_spar_from_cavity(0.004, 0.004, gap_min=0.002) == 0.0


@pytest.mark.parametrize("d", [0.003, 0.00667, 0.010, 0.020, 0.100])
def test_cavity_bound_with_defaults_is_half_the_clamp_allowance(d: float) -> None:
    """SG-039: with the default floor the gate equals the clamp's allowance (R10.7) / 2."""
    from implicit_geometry_kernel.section.assembly.cap_feasibility import (
        MIN_CAP_VOID_ABS_M,
        MIN_CAP_VOID_FRAC,
    )

    t_s = 0.004
    allowance = d - max(MIN_CAP_VOID_ABS_M, MIN_CAP_VOID_FRAC * d)
    assert max_spar_from_cavity(d + t_s, t_s) == pytest.approx(0.5 * allowance, abs=1e-15)


@pytest.mark.parametrize("h", [0.0, -1.0, math.nan])
def test_cavity_bound_is_unbounded_without_a_cavity_height(h: float) -> None:
    assert max_spar_from_cavity(h, 0.004) == math.inf


def test_a_thin_cavity_steps_the_layout_to_its_capless_sibling() -> None:
    layout = resolve_system_type("2D-F-LC-0R")
    out = apply_spar_cap_runout(
        _Resolution(layout, (0.2, 0.5)),
        h_local=0.010,
        gap_min=0.002,
        skin_thickness_floor=0.004,
        spar_min_gauge=0.005,
    )
    assert not out.layout.has_spar_caps
    assert out.web_x_positions_m == (0.2, 0.5)


def test_a_deep_cavity_keeps_the_caps() -> None:
    layout = resolve_system_type("2D-F-LC-0R")
    res = _Resolution(layout, (0.2, 0.5))
    out = apply_spar_cap_runout(
        res, h_local=0.20, gap_min=0.002, skin_thickness_floor=0.004, spar_min_gauge=0.005
    )
    assert out is res


@pytest.mark.parametrize("key", ["2D-F-LC-0R", "1C-F-LC-0R"])
def test_placement_check_passes_on_a_thick_airfoil(key: str) -> None:
    layout = resolve_system_type(key)
    xs = (0.2, 0.5) if key.startswith("2") else (0.3,)
    res = _Resolution(layout, xs)
    out = apply_spar_cap_runout(
        res,
        h_local=None,
        gap_min=0.002,
        skin_thickness_floor=0.004,
        spar_min_gauge=0.005,
        airfoil=_naca("0024", chord=2.0),
    )
    assert out is res


def test_placement_check_runs_out_caps_that_cannot_fit() -> None:
    layout = resolve_system_type("2D-F-LC-0R")
    out = apply_spar_cap_runout(
        _Resolution(layout, (0.2, 0.5)),
        h_local=None,
        gap_min=0.002,
        skin_thickness_floor=0.004,
        spar_min_gauge=0.05,
        airfoil=_naca("0012", chord=0.5),
    )
    assert not out.layout.has_spar_caps


def test_capless_layouts_zero_gauge_and_none_pass_through() -> None:
    capless = _Resolution(resolve_system_type("2A-F-LC-0R"), (0.2, 0.5))
    assert apply_spar_cap_runout(
        capless, h_local=0.001, gap_min=0.002, skin_thickness_floor=0.004, spar_min_gauge=0.005
    ) is capless
    capped = _Resolution(resolve_system_type("2D-F-LC-0R"), (0.2, 0.5))
    assert apply_spar_cap_runout(
        capped, h_local=0.001, gap_min=0.002, skin_thickness_floor=0.004, spar_min_gauge=0.0
    ) is capped
    assert apply_spar_cap_runout(
        None, h_local=0.001, gap_min=0.002, skin_thickness_floor=0.004, spar_min_gauge=0.005
    ) is None


def _fixed_band_layout(**rule_kw):
    import dataclasses

    from implicit_geometry_kernel.section.assembly.from_layout import as_section_layout
    from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule

    base = as_section_layout(resolve_system_type("2B-F-LC-0R"))
    return dataclasses.replace(base, cap_placement=CapPlacementRule.fixed_band(**rule_kw))


def test_fixed_band_gate_tests_the_band_the_build_clamps() -> None:
    """SG-041: the band is centred on the rule's anchor, not on mid-chord."""
    from implicit_geometry_kernel.section.layout.section_axes import max_thickness_chord_x
    from implicit_geometry_kernel.section.layout.spar_cap_runout import _fixed_bands

    af = _naca("2412")
    rule = _fixed_band_layout(anchor="max_thickness", half_width=0.1).cap_placement
    bands, anchor = _fixed_bands(rule, af, 1.0, None)
    x_mt = max_thickness_chord_x(af)
    assert anchor == "band_center"
    assert bands == [pytest.approx((x_mt - 0.1, x_mt + 0.1))]
    pitching = _fixed_band_layout(anchor="pitching", half_width=0.1).cap_placement
    with pytest.warns(DeprecationWarning, match="pitch_fraction"):
        legacy, _ = _fixed_bands(pitching, af, 1.0, None)
    assert legacy == [pytest.approx((0.4, 0.6))]
    bands, _ = _fixed_bands(pitching, af, 1.0, 0.35)
    assert 0.5 * (bands[0][0] + bands[0][1]) == pytest.approx(0.35, abs=1e-3)


def test_pad_gate_uses_the_given_half_width() -> None:
    """SG-041: a wide pad that reaches past the cavity runs out; the default 0.04c pad fits."""
    layout = resolve_system_type("1C-F-LC-0R")
    res = _Resolution(layout, (0.85,))
    kw = dict(h_local=None, gap_min=0.002, skin_thickness_floor=0.004, spar_min_gauge=0.005,
              airfoil=_naca("2412"))
    assert apply_spar_cap_runout(res, **kw) is res
    assert not apply_spar_cap_runout(res, pad_half_width=0.3, **kw).layout.has_spar_caps
