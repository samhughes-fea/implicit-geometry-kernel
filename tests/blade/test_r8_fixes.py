"""Regressions for the loft defects found while deriving R8 (BG-058 to BG-060)."""

from __future__ import annotations

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.loft.span_place import offset_rings
from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout


def _blade(af, **kw):
    stations = [Station(z=0.0, chord=1.0, airfoil=af), Station(z=2.0, chord=1.0, airfoil=af)]
    kw.setdefault("pivot_fraction", 0.25)
    return BladeGeometry(stations, kw.pop("layout", "2D-F-LC-0R"), kw.pop("skin_thickness", 0.004), **kw)


@pytest.mark.parametrize("skin", [0.004, 0.001])
def test_thin_walled_cut_keeps_metre_web_stations(skin):
    """BG-058: web_stations_m is scaled with the chord when the cut is built k times larger."""
    layout = SectionLayout(
        web_chord_fracs=(0.3, 0.6),
        orientation="vertical",
        cap_placement=CapPlacementRule.none(),
        web_stations_m=(0.3, 0.6),
    )
    cut = _blade(Airfoil.naca4("0018"), layout=layout, skin_thickness=skin, shear_web_thickness=0.005).section_at(1.0)
    xs = sorted(float(np.mean(w.points[:, 0])) + 0.25 for w in cut.walls if w.kind == "web")
    assert xs == pytest.approx([0.3, 0.6], abs=1e-9)


def test_scale_layout_metres_scales_every_metre_field():
    """BG-058: stations, axes and a fixed cap half-width all scale; fractions do not."""
    from implicit_geometry_kernel.blade.loft.loft_cut import _scale_layout_metres

    layout = SectionLayout(
        web_chord_fracs=(0.3,),
        orientation="vertical",
        cap_placement=CapPlacementRule.fixed_band(anchor="max_thickness", half_width=0.1),
        web_stations_m=(0.3,),
        web_axes_m=(((0.3, 0.05), (0.31, -0.04)),),
    )
    out = _scale_layout_metres(layout, 2.5)
    assert out.web_chord_fracs == (0.3,)
    assert out.web_stations_m == pytest.approx((0.75,))
    assert np.allclose(np.asarray(out.web_axes_m), 2.5 * np.asarray(layout.web_axes_m))
    assert out.cap_placement.half_width == pytest.approx(0.25)
    assert _scale_layout_metres(layout, 1.0) is layout


def test_outer_cap_edge_moves_inward_on_a_pressure_surface_above_the_chord():
    """BG-059: the inward normal comes from the side, not from a point at mid-chord."""
    geo = _blade(
        Airfoil.naca4("4406"),
        skin_thickness=0.01,
        spar_cap_thickness=0.04,
        shear_web_thickness=0.005,
        cap_edge_reference="outer",
        cap_windows=lambda z: (0.2, 0.5),
    )
    _, _, _, _, (x, _, lower) = geo.interpolate(1.0)
    f = 0.2
    assert float(np.interp(f, x, lower)) > 0.0  # the pressure surface is above the chord here
    slope = float(np.interp(f, x, np.gradient(lower, x)))
    inward_xi = -slope / np.hypot(1.0, slope)
    expected = f + inward_xi * (0.01 + 0.02) / 1.0
    assert geo._outer_fraction_to_midplane(1.0, f, "pressure") == pytest.approx(expected, abs=1e-4)


@pytest.mark.parametrize("code", ["0012", "4412", "4406", "6409"])
@pytest.mark.parametrize("reverse", [False, True])
def test_offset_rings_moves_every_vertex_inward(code, reverse):
    """BG-060: one orientation per ring from its signed area, either winding."""
    foil = Airfoil.naca4(code).to_sdf(chord=1.0)
    xy = np.asarray(foil.vertices, float)
    if np.allclose(xy[0], xy[-1]):
        xy = xy[:-1]
    if reverse:
        xy = xy[::-1]
    ring = np.column_stack([xy, np.zeros(len(xy))])
    ring = np.vstack([ring, ring[:1]])[None]
    moved = offset_rings(ring, 0.002)[0, :-1, :2]
    a2 = np.sum(xy[:, 0] * np.roll(xy[:, 1], -1) - np.roll(xy[:, 0], -1) * xy[:, 1])
    tang = np.roll(xy, -1, 0) - np.roll(xy, 1, 0)
    inward = np.column_stack([-tang[:, 1], tang[:, 0]]) * np.sign(a2)
    assert np.all(np.sum(inward * (moved - xy), axis=1) > 0.0)
