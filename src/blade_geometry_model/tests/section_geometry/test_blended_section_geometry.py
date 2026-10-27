"""Unit tests for BlendedSectionGeometry.

Covers:
- touch_points_S returns all four NESW keys with correct shapes
- polygon_S returns correctly shaped geometry and counts
- x_scale_factors scales only the X_S (chord-normal) axis
- station_indices selects a subset of stations
- polygon_S pads to a requested max_pts larger than native
"""
from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.root_blend_model import (
    ExactAirfoilGeometry,
    RootBlendConfig,
    RootBlendGridSpec,
    SpanwiseRegions,
    QuadTreeSpec,
    build_root_blend_stations,
    exact_airfoil_geometry_from_naca,
)
from blade_geometry_model.section_models import BlendedSectionGeometry


# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

_N = 8  # number of stations


def _make_stations():
    """Build a minimal RootBlendStations for testing."""
    regions = SpanwiseRegions(
        eta_cyl_start=0.125,
        eta_cyl_end=0.20,
        eta_tr_start=0.20,
        eta_tr_end=0.40,
        eta_air_start=0.40,
        eta_air_end=1.0,
    )
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 151)))
    airfoil = exact_airfoil_geometry_from_naca(x, thickness_ratio=0.35)
    cfg = RootBlendConfig(
        regions=regions,
        radius_over_chord_inboard=0.35,
        airfoil=airfoil,
        chord=np.full(_N, 2.0),
        pitch_axis_over_chord=np.linspace(0.50, 0.333, _N),
        quadtree=QuadTreeSpec(
            x_min=-0.8, x_max=0.8, y_min=-0.5, y_max=0.5,
            max_depth=4, min_cell_size=0.02, refine_band=0.05, max_cells=4000,
        ),
        grid_spec=RootBlendGridSpec(nx=181, ny=181),
    )

    class _Dims:
        pass

    d = _Dims()
    d.chord = np.full(_N, 2.0)
    d.twist = np.linspace(0.0, 5.0, _N)
    d.outer_loop_geometry_S = np.zeros((_N, 256, 2))
    return build_root_blend_stations(d, cfg)


@pytest.fixture(scope="module")
def stations():
    return _make_stations()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestTouchPointsS:
    def test_returns_all_keys(self, stations):
        geom = BlendedSectionGeometry(stations)
        tp = geom.touch_points_S(_N)
        assert set(tp.keys()) == {"north", "south", "east", "west"}

    def test_shape(self, stations):
        geom = BlendedSectionGeometry(stations)
        tp = geom.touch_points_S(_N)
        for key, arr in tp.items():
            assert arr.shape == (_N, 2), f"key={key}: expected ({_N}, 2), got {arr.shape}"


class TestPolygonS:
    def test_shape_at_native_max_pts(self, stations):
        geom = BlendedSectionGeometry(stations)
        native_max = stations.outer_loop_geometry_S.shape[1]
        poly, counts = geom.polygon_S(_N, native_max)
        assert poly.shape == (_N, native_max, 2)
        assert counts.shape == (_N,)
        assert np.all(counts <= native_max)

    def test_padding_to_larger_max_pts(self, stations):
        geom = BlendedSectionGeometry(stations)
        native_max = stations.outer_loop_geometry_S.shape[1]
        big = native_max + 50
        poly, counts = geom.polygon_S(_N, big)
        assert poly.shape == (_N, big, 2)
        # Padded tail should be zero
        for i in range(_N):
            n = int(counts[i])
            assert np.all(poly[i, n:, :] == 0.0), f"Station {i}: non-zero padding beyond count"

    def test_truncation_to_smaller_max_pts(self, stations):
        geom = BlendedSectionGeometry(stations)
        small = 32
        poly, counts = geom.polygon_S(_N, small)
        assert poly.shape == (_N, small, 2)
        assert np.all(counts <= small)


class TestXScaleFactors:
    def test_north_south_x_scaled(self, stations):
        """X_S (chord-normal, col 0) on north/south should be multiplied by the scale."""
        base = BlendedSectionGeometry(stations)
        scaled = BlendedSectionGeometry(stations, x_scale_factors=np.full(_N, 2.0))
        tp_base = base.touch_points_S(_N)
        tp_scaled = scaled.touch_points_S(_N)
        np.testing.assert_allclose(
            tp_scaled["north"][:, 0], 2.0 * tp_base["north"][:, 0], atol=1e-12
        )
        np.testing.assert_allclose(
            tp_scaled["south"][:, 0], 2.0 * tp_base["south"][:, 0], atol=1e-12
        )

    def test_east_west_y_not_scaled(self, stations):
        """Y_S (chordwise, col 1) on east/west must be unchanged by x_scale_factors."""
        base = BlendedSectionGeometry(stations)
        scaled = BlendedSectionGeometry(stations, x_scale_factors=np.full(_N, 3.0))
        tp_base = base.touch_points_S(_N)
        tp_scaled = scaled.touch_points_S(_N)
        np.testing.assert_allclose(tp_scaled["east"], tp_base["east"], atol=1e-12)
        np.testing.assert_allclose(tp_scaled["west"], tp_base["west"], atol=1e-12)

    def test_polygon_x_column_scaled(self, stations):
        """polygon_S col 0 (X_S) within active points should be multiplied by scale."""
        scale = 1.5
        base = BlendedSectionGeometry(stations)
        scaled = BlendedSectionGeometry(stations, x_scale_factors=np.full(_N, scale))
        native_max = stations.outer_loop_geometry_S.shape[1]
        poly_b, counts = base.polygon_S(_N, native_max)
        poly_s, _ = scaled.polygon_S(_N, native_max)
        for i in range(_N):
            n = int(counts[i])
            np.testing.assert_allclose(
                poly_s[i, :n, 0], scale * poly_b[i, :n, 0], atol=1e-12,
                err_msg=f"Station {i}: X_S column not scaled correctly",
            )
            np.testing.assert_allclose(
                poly_s[i, :n, 1], poly_b[i, :n, 1], atol=1e-12,
                err_msg=f"Station {i}: Y_S column should not be scaled",
            )


class TestStationIndices:
    def test_subset_leading_dimension(self, stations):
        idx = np.array([2, 4])
        geom = BlendedSectionGeometry(stations, station_indices=idx)
        tp = geom.touch_points_S(len(idx))
        for key, arr in tp.items():
            assert arr.shape[0] == 2, f"key={key}: expected 2 rows, got {arr.shape[0]}"

    def test_subset_matches_full_at_same_indices(self, stations):
        idx = np.array([1, 5])
        full = BlendedSectionGeometry(stations)
        subset = BlendedSectionGeometry(stations, station_indices=idx)
        tp_full = full.touch_points_S(_N)
        tp_sub = subset.touch_points_S(len(idx))
        for key in ("north", "south", "east", "west"):
            np.testing.assert_allclose(
                tp_sub[key], tp_full[key][idx], atol=1e-12,
                err_msg=f"Subset mismatch for key={key}",
            )
