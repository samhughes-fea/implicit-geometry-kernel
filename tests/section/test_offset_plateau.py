"""Offset-branch plateaus are collapsed; true reversals stay infeasible."""

import numpy as np

from implicit_geometry_kernel.section.assembly.clips import _airfoil_offset_branch, _strictly_x_ascending
from implicit_geometry_kernel.sdf.airfoil import AirfoilSDF


def test_collapse_duplicate_x():
    br = np.array([[0.0, 0.0], [1.0, 0.1], [1.0, 0.1], [2.0, 0.0]])
    out = _strictly_x_ascending(br)
    assert np.all(np.diff(out[:, 0]) > 0.0)
    assert len(out) == 3


def test_naca_offset_branch_is_graph():
    verts = np.vstack([
        np.column_stack([np.linspace(1.0, 0.0, 40), 0.12 * np.sin(np.linspace(0, np.pi, 40))]),
        np.column_stack([np.linspace(0.0, 1.0, 40)[1:], -0.12 * np.sin(np.linspace(0, np.pi, 40))[1:]]),
    ])
    foil = AirfoilSDF(verts, chord=1.0)
    upper = _airfoil_offset_branch(foil, depth=0.01, surface="upper")
    assert np.all(np.diff(upper[:, 0]) > 0.0)


def test_true_reversal_still_raises():
    br = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.1]])
    out = _strictly_x_ascending(br)
    assert np.any(np.diff(out[:, 0]) < 0.0)
