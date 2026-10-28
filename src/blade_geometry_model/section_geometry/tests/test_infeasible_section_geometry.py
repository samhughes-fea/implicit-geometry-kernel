"""Infeasible sections raise a type the sizing loop can catch.

A thickness search will inevitably propose sections that cannot be built — a cap
deeper than the local cavity, a skin that closes the section, a web with no
inner-skin intercept. Those must be *scored* as failed trials so the search backs
off, not raised through the station worker.

This used to be decided by substring-matching the exception message in a BSP
optimize module, which meant any newly-worded or newly-added infeasibility became
a crash. These tests pin the type instead. BSP sizing absorption is covered in
``blade_structure_precompute.module.section_selection`` (station sizing workflow);
this package only owns the geometry oracles.
"""

from __future__ import annotations

import pytest

from blade_geometry_model.section_contract import (
    InfeasibleSectionGeometry,
)


def test_is_a_valueerror_subclass() -> None:
    """Pre-existing ``except ValueError`` handlers must keep working."""
    assert issubclass(InfeasibleSectionGeometry, ValueError)


def test_offset_ring_deeper_than_the_section_is_infeasible() -> None:
    """The direct oracle: no locus exists at a depth past the inscribed limit."""
    from blade_geometry_model.section_geometry.constructive.offset_locus import (
        airfoil_offset_ring,
    )
    from blade_geometry_model.section_geometry.sdf.airfoil import (
        AirfoilSDF,
    )

    # chord 0.112 m has max inscribed depth ~7.84 mm.
    af = AirfoilSDF.from_naca_series(
        6, 63.0, 4.0, 15.0, n_points=200, chord=0.112, closed_te=True
    )
    airfoil_offset_ring(af, 0.003)  # feasible, must not raise

    with pytest.raises(InfeasibleSectionGeometry):
        airfoil_offset_ring(af, 0.009)


def test_bad_depth_argument_is_not_infeasibility() -> None:
    """A malformed request is a caller bug, not an unbuildable section.

    Keeping these distinct is the point: the sizing loop must not swallow real
    programming errors as failed trials.
    """
    from blade_geometry_model.section_geometry.constructive.offset_locus import (
        airfoil_offset_ring,
    )
    from blade_geometry_model.section_geometry.sdf.airfoil import (
        AirfoilSDF,
    )

    af = AirfoilSDF.from_naca("0012", n_points=100, chord=1.0)
    with pytest.raises(ValueError) as exc:
        airfoil_offset_ring(af, -0.001)
    assert not isinstance(exc.value, InfeasibleSectionGeometry)


def test_web_anchors_and_skin_intercepts_are_infeasibility() -> None:
    """The two cases the old message allowlist covered, now matched by type."""
    from blade_geometry_model.section_geometry.sections.web_layout import (
        inner_y_at_x,
    )
    from blade_geometry_model.section_geometry.sdf.airfoil import (
        AirfoilSDF,
    )

    af = AirfoilSDF.from_naca("0012", n_points=100, chord=1.0)
    # A skin thicker than the local half-thickness leaves no inner cavity.
    with pytest.raises(InfeasibleSectionGeometry):
        inner_y_at_x(af, 0.4, 0.95)
