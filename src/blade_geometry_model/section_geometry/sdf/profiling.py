"""Compatibility re-export — prefer ``section_geometry.profiling``."""

from blade_geometry_model.section_geometry.profiling import *  # noqa: F403
from blade_geometry_model.section_geometry.profiling import (
    SdfProfilingStats,
    get_stats,
    profiling_enabled,
    record_grid_eval,
    record_sdf_polygon,
    reset_stats,
    timed_block,
    timed_stage,
)

__all__ = [
    "SdfProfilingStats",
    "get_stats",
    "profiling_enabled",
    "record_grid_eval",
    "record_sdf_polygon",
    "reset_stats",
    "timed_block",
    "timed_stage",
]
