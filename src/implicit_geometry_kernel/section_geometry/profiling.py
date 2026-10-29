"""Optional section_geometry profiling (``SECTION_GEOMETRY_PROFILE=1``).

Cross-cutting timing helpers shared by api / assemble / constructive / topology.
Kept out of ``sdf/`` so G-row homes do not skip-import the SDF kernel for timers.
"""

from __future__ import annotations

import os
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any


def profiling_enabled() -> bool:
    return os.getenv("SECTION_GEOMETRY_PROFILE", "0").lower() in {"1", "true", "yes"}


@dataclass
class SdfProfilingStats:
    """Thread-local-friendly counters for section_geometry hot paths."""

    grid_eval_count: int = 0
    sdf_polygon_count: int = 0
    sdf_polygon_points: int = 0
    sdf_polygon_s: float = 0.0
    station_build_s: float = 0.0
    properties_s: float = 0.0
    plot_s: float = 0.0
    phi_eval_count: int = 0
    quadtree_leaf_count: int = 0
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "grid_eval_count": int(self.grid_eval_count),
            "sdf_polygon_count": int(self.sdf_polygon_count),
            "sdf_polygon_points": int(self.sdf_polygon_points),
            "sdf_polygon_s": float(self.sdf_polygon_s),
            "phi_eval_count": int(self.phi_eval_count),
            "quadtree_leaf_count": int(self.quadtree_leaf_count),
            "station_build_s": float(self.station_build_s),
            "properties_s": float(self.properties_s),
            "plot_s": float(self.plot_s),
            **{str(k): v for k, v in self.extra.items()},
        }


_local = threading.local()


def get_stats() -> SdfProfilingStats | None:
    if not profiling_enabled():
        return None
    stats = getattr(_local, "stats", None)
    if stats is None:
        stats = SdfProfilingStats()
        _local.stats = stats
    return stats


def reset_stats() -> SdfProfilingStats | None:
    if not profiling_enabled():
        return None
    stats = SdfProfilingStats()
    _local.stats = stats
    return stats


@contextmanager
def timed_block(name: str):
    stats = get_stats()
    if stats is None:
        yield
        return
    t0 = time.perf_counter()
    try:
        yield
    finally:
        dt = time.perf_counter() - t0
        if name == "station_build":
            stats.station_build_s += dt
        elif name == "properties":
            stats.properties_s += dt
        elif name == "plot":
            stats.plot_s += dt
        else:
            stats.extra[name] = float(stats.extra.get(name, 0.0)) + dt


@contextmanager
def timed_stage(name: str, sink: dict[str, float] | None = None):
    """Time a named stage into ``sink`` and, when enabled, into profiling extras.

    ``sink`` is filled even when ``SECTION_GEOMETRY_PROFILE`` is off so callers
    can decide whether to attach timings. The reserved ``properties`` /
    ``station_build`` / ``plot`` buckets are left alone — ``props_compute``
    already records ``properties_s`` via :func:`timed_block`.
    """
    t0 = time.perf_counter()
    try:
        yield
    finally:
        dt = time.perf_counter() - t0
        if sink is not None:
            sink[name] = float(sink.get(name, 0.0)) + dt
        stats = get_stats()
        if stats is not None and name not in {"properties", "station_build", "plot"}:
            stats.extra[name] = float(stats.extra.get(name, 0.0)) + dt


def record_grid_eval(n_points: int = 0) -> None:
    stats = get_stats()
    if stats is not None:
        stats.grid_eval_count += 1
        if n_points:
            stats.phi_eval_count += int(n_points)


def record_sdf_polygon(n_points: int, elapsed_s: float = 0.0) -> None:
    stats = get_stats()
    if stats is not None:
        stats.sdf_polygon_count += 1
        stats.sdf_polygon_points += int(n_points)
        if elapsed_s:
            stats.sdf_polygon_s += float(elapsed_s)
