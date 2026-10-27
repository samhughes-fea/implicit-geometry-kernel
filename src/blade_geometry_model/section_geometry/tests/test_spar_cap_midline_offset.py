"""The spar-cap midline must lie on its own solid's mid-surface.

``SparCap`` builds the solid as ``Slab(Offset(poly, -t_skin), cap_height)`` — the
band between the true *normal* offsets at depths ``t_skin`` and
``t_skin + cap_height`` — so the mid-surface is the normal offset at
``t_skin + cap_height/2``. The midline was previously a *vertical* displacement
of the outer surface, which lands at normal depth ``d*cos(theta)`` and therefore
ran shallow wherever the surface is steep; in the worst measured case it sat
inside the skin band rather than in the cap, with 86 % depth error.

Nothing pinned that before, which is why it went unnoticed. These tests pin it.
Segment midpoints are sampled as well as vertices: sampling vertices alone
passes even for a polyline that chords across the locus' corners, which is the
part that is easy to get wrong.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.sdf.primitives import sdf_polygon

# Chord-relative thicknesses, so no case is asked for a cap that cannot fit.
_T_SKIN_FRAC = 0.003
_T_CAP_FRAC = 0.012

AIRFOILS = {
    "0012@1.0": lambda: AirfoilSDF.from_naca("0012", n_points=200, chord=1.0),
    "2412@1.0": lambda: AirfoilSDF.from_naca("2412", n_points=200, chord=1.0),
    "6s@1.656": lambda: AirfoilSDF.from_naca_series(
        6, 63.0, 4.0, 15.0, n_points=200, chord=1.656, closed_te=True
    ),
    "6s@0.241": lambda: AirfoilSDF.from_naca_series(
        6, 63.0, 4.0, 15.0, n_points=200, chord=0.241, closed_te=True
    ),
}
LAYOUTS = ["2D-F", "2D-CN", "2C-CN", "3D-F", "1C-CN", "5D-CN"]


def _caps(airfoil_key: str, layout_key: str, twist_rad: float = 0.0):
    af = AIRFOILS[airfoil_key]()
    chord = float(getattr(af, "chord", 1.0))
    section = build_section_view(
        af,
        resolve_system_type(layout_key),
        twist_angle_rad=twist_rad,
        skin_thickness=_T_SKIN_FRAC * chord,
        web_thickness=_T_SKIN_FRAC * chord,
        cap_height=_T_CAP_FRAC * chord,
    )
    caps = getattr(section, "spar_cap_components_unrotated", {}) or {}
    return af, caps


def _sample(arc: np.ndarray) -> np.ndarray:
    """Vertices plus segment midpoints."""
    return np.vstack([arc, 0.5 * (arc[:-1] + arc[1:])])


@pytest.mark.parametrize("airfoil_key", list(AIRFOILS))
@pytest.mark.parametrize("layout_key", LAYOUTS)
def test_cap_midline_is_at_normal_offset_depth(airfoil_key: str, layout_key: str) -> None:
    """phi_airfoil along the midline must equal -(t_skin + cap_h/2)."""
    af, caps = _caps(airfoil_key, layout_key)
    verts = np.asarray(af.vertices, dtype=float)
    assert caps, f"{layout_key} produced no spar caps"

    for label, comp in caps.items():
        arc = np.asarray(comp.midline_polyline(n=40), dtype=float)
        if arc.shape[0] < 2:
            continue
        depth = float(comp.skin_thickness) + float(comp.cap_height) / 2.0
        pts = _sample(arc)
        phi = np.asarray(sdf_polygon(pts[:, 0], pts[:, 1], verts), dtype=float)
        err = float(np.max(np.abs(phi + depth))) / depth
        assert err < 5e-3, (
            f"{airfoil_key} {layout_key} {label}: midline is {err:.3%} off the "
            f"normal-offset depth {depth:g} m — it is not on the cap solid's "
            "mid-surface"
        )


@pytest.mark.parametrize("airfoil_key", list(AIRFOILS))
@pytest.mark.parametrize("layout_key", LAYOUTS)
def test_cap_midline_lies_inside_the_cap_slab(airfoil_key: str, layout_key: str) -> None:
    """Clip-independent form of the same claim: inside the unclipped slab."""
    from blade_geometry_model.section_geometry.sections.clip_algebra import (
        airfoil_expr,
    )
    from blade_geometry_model.section_geometry.recipes.sdf_ir.base import (
        expr_callable,
    )
    from blade_geometry_model.section_geometry.sdf.ir import Offset, Slab

    af, caps = _caps(airfoil_key, layout_key)
    poly = airfoil_expr(af)

    for label, comp in caps.items():
        arc = np.asarray(comp.midline_polyline(n=40), dtype=float)
        if arc.shape[0] < 2:
            continue
        cap_h = float(comp.cap_height)
        slab = expr_callable(
            Slab(Offset(poly, amount=-float(comp.skin_thickness)), depth=cap_h)
        )
        pts = _sample(arc)
        phi = np.asarray(slab(pts[:, 0], pts[:, 1]), dtype=float)
        # Well inside, not merely within: the mid-surface should sit ~cap_h/2
        # from either face.
        assert float(np.max(phi)) <= -0.4 * cap_h, (
            f"{airfoil_key} {layout_key} {label}: midline point at phi_slab="
            f"{float(np.max(phi)):.3e} is not comfortably inside the cap band "
            f"(cap_height={cap_h:g} m)"
        )


@pytest.mark.parametrize("airfoil_key", list(AIRFOILS))
@pytest.mark.parametrize("layout_key", LAYOUTS)
def test_cap_midline_satisfies_strip_contract(airfoil_key: str, layout_key: str) -> None:
    """Ascending x, non-degenerate segments, and the documented endpoints."""
    from blade_geometry_model.section_mesh.handoff.shell_inputs import (
        _MIN_SEGMENT_LENGTH_M,
    )

    _af, caps = _caps(airfoil_key, layout_key)

    for label, comp in caps.items():
        arc = np.asarray(comp.midline_polyline(n=40), dtype=float)
        if arc.shape[0] < 2:
            continue
        dx = np.diff(arc[:, 0])
        assert np.all(dx > 0.0), (
            f"{airfoil_key} {layout_key} {label}: midline is not strictly "
            "ascending in chord-frame x, violating the ShellMidlineStrip contract"
        )
        seg = np.hypot(np.diff(arc[:, 0]), np.diff(arc[:, 1]))
        assert float(seg.min()) > _MIN_SEGMENT_LENGTH_M, (
            f"{airfoil_key} {layout_key} {label}: shortest segment "
            f"{float(seg.min()):.3e} m would be rejected by _validate_strips"
        )
        if getattr(comp, "_strip_clip", None) is None:
            span = float(comp.x_end) - float(comp.x_start)
            assert abs(arc[0, 0] - float(comp.x_start)) < 1e-9 * max(span, 1e-9)
            assert abs(arc[-1, 0] - float(comp.x_end)) < 1e-9 * max(span, 1e-9)
        else:
            phi_end = np.asarray(
                comp._strip_clip(arc[[0, -1], 0], arc[[0, -1], 1]), dtype=float
            )
            scale = max(abs(float(comp.x_end) - float(comp.x_start)), 1e-9)
            assert np.all(np.abs(phi_end) < 1e-6 * scale), (
                f"{airfoil_key} {layout_key} {label}: clipped endpoints are not "
                f"on the strip boundary (phi={phi_end})"
            )


def test_vertical_offset_would_fail_the_depth_gate() -> None:
    """Guard the gate itself: the superseded construction must not pass it.

    Without this, a regression back to a vertical offset could go unnoticed if
    the tolerance above were ever loosened.
    """
    af = AIRFOILS["6s@0.241"]()
    chord = float(getattr(af, "chord", 1.0))
    verts = np.asarray(af.vertices, dtype=float)
    t_skin, cap_h = _T_SKIN_FRAC * chord, _T_CAP_FRAC * chord
    depth = t_skin + cap_h / 2.0

    seg = np.asarray(af.upper_surface(), dtype=float)
    seg = seg[np.argsort(seg[:, 0])]
    xs = np.linspace(0.25 * chord, 0.55 * chord, 40)
    ys = np.interp(xs, seg[:, 0], seg[:, 1]) - depth  # the old vertical offset
    phi = np.asarray(sdf_polygon(xs, ys, verts), dtype=float)
    err = float(np.max(np.abs(phi + depth))) / depth
    assert err > 5e-3, (
        "the vertical-offset construction now passes the depth gate; the gate "
        "no longer discriminates and its tolerance needs revisiting"
    )
