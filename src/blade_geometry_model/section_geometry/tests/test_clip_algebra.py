"""Shared clip algebra IR + Shapely round-trip smoke tests."""

from __future__ import annotations

from shapely.geometry import box

from blade_geometry_model.section_geometry.sdf.ir import Intersect
from blade_geometry_model.section_geometry.sections.clip_algebra import (
    ChordBand,
    apply_shapely,
    to_ir,
)


def test_chord_band_ir_and_shapely_reduces_x_extent():
    clip = ChordBand(x_start=0.2, x_end=0.6)
    ir = to_ir(clip)
    assert isinstance(ir, Intersect)

    poly = box(0.0, -0.1, 1.0, 0.1)
    clipped = apply_shapely(clip, poly)
    assert not clipped.is_empty
    minx, _miny, maxx, _maxy = clipped.bounds
    assert minx >= 0.2 - 1e-9
    assert maxx <= 0.6 + 1e-9
    assert (maxx - minx) < (1.0 - 0.0) - 1e-9