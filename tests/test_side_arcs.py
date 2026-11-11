"""Each side arc runs from its own TE corner to the LE without passing the other corner."""

import numpy as np

from implicit_geometry_kernel.blade.topology.map_walk import _side_arcs


def test_pressure_arc_does_not_walk_through_the_suction_corner():
    # the LE sits nearer the pressure corner the short way round, which passes
    # the suction corner (a blunt TE whose corners are a few vertices apart)
    n = 20
    t = 2.0 * np.pi * np.arange(n) / n
    pts = np.column_stack([np.cos(t), np.sin(t), np.zeros(n)])
    su, pr = _side_arcs(pts, None, i_le=14, i_su=1, i_pr=3)
    assert np.array_equal(su[0], pts[1]) and np.array_equal(pr[0], pts[3])
    assert np.array_equal(su[-1], pts[14]) and np.array_equal(pr[-1], pts[14])
    assert not any(np.array_equal(p, pts[3]) for p in su)
    assert not any(np.array_equal(p, pts[1]) for p in pr)
    # every loop vertex once, except the TE wrap between the corners; the LE twice
    assert len(su) + len(pr) == n - 1 + 1
