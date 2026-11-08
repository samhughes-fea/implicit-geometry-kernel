"""Shared clip algebra IR + polygon round-trip smoke tests."""

from __future__ import annotations

from implicit_geometry_kernel.sdf.nodes import Intersect
from implicit_geometry_kernel.section.assembly.clip_algebra import (
    ChordBand,
    apply_clip,
    to_ir,
)
from implicit_geometry_kernel.section.boundary import envelope


def test_chord_band_ir_and_clip_reduces_x_extent():
    clip = ChordBand(x_start=0.2, x_end=0.6)
    ir = to_ir(clip)
    assert isinstance(ir, Intersect)

    poly = envelope(0.0, -0.1, 1.0, 0.1)
    clipped = apply_clip(clip, poly)
    assert not clipped.is_empty
    minx, _miny, maxx, _maxy = clipped.bounds
    assert minx >= 0.2 - 1e-9
    assert maxx <= 0.6 + 1e-9
    assert (maxx - minx) < (1.0 - 0.0) - 1e-9