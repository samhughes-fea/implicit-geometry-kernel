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


def test_cavity_bound_is_half_the_room_left_by_skin_and_gap() -> None:
    assert max_spar_from_cavity(0.10, 0.004, gap_min=0.002) == pytest.approx(0.047)
    assert max_spar_from_cavity(0.004, 0.004, gap_min=0.002) == 0.0


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
