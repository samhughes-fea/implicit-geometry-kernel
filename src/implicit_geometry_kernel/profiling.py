"""Optional section_geometry profiling (``SECTION_GEOMETRY_PROFILE=1``).

Cross-cutting timing helpers shared by api / assemble / constructive / topology.
Kept out of ``sdf/`` so section packages do not skip-import the SDF kernel for timers.

Wall-clock stages report on the shared ``blade.stage`` logger. A caller attaches
:class:`StageCollector` and folds the reports into one timings object. A stage
that is not entered is absent. A stage that raises reports the seconds up to
the error and the error text, then the exception propagates. The level sits
above INFO, so silencing INFO does not drop a stage report. Times are seconds
from :func:`time.perf_counter`.
"""

from __future__ import annotations

import logging
import os
import threading
import time
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from contextvars import ContextVar
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
    already records ``properties_s`` via :func:`timed_block`. The same seconds
    are reported on the shared ``blade.stage`` log.
    """
    t0 = time.perf_counter()
    error: str | None = None
    try:
        yield
    except Exception as exc:
        error = str(exc)
        raise
    finally:
        dt = time.perf_counter() - t0
        if sink is not None:
            sink[name] = float(sink.get(name, 0.0)) + dt
        stats = get_stats()
        if stats is not None and name not in {"properties", "station_build", "plot"}:
            stats.extra[name] = float(stats.extra.get(name, 0.0)) + dt
        report_stage("implicit_geometry_kernel", name, dt, error, t0)


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


STAGE = 25
LOGGER_NAME = "blade.stage"
logging.addLevelName(STAGE, "STAGE")

_depth: ContextVar[int] = ContextVar("blade_stage_depth", default=0)


@dataclass(frozen=True)
class StageReport:
    """One finished stage."""

    package: str
    stage: str
    seconds: float
    error: str | None
    depth: int
    started: float


def report_stage(package: str, name: str, seconds: float, error: str | None,
                 started: float, depth: int | None = None) -> None:
    """Emit one stage report. ``depth`` defaults to the current nesting level."""
    level = _depth.get() if depth is None else depth
    logging.getLogger(LOGGER_NAME).log(
        STAGE,
        "%s %s %.6fs%s",
        package,
        name,
        seconds,
        "" if error is None else f" error={error}",
        extra={
            "stage_package": package,
            "stage_name": name,
            "stage_seconds": float(seconds),
            "stage_error": error,
            "stage_depth": int(level),
            "stage_started": float(started),
        },
    )


@contextmanager
def stage(package: str, name: str) -> Iterator[None]:
    """Time ``name`` and report it when the block ends."""
    depth = _depth.get()
    token = _depth.set(depth + 1)
    started = time.perf_counter()
    try:
        yield
    except Exception as exc:
        report_stage(package, name, time.perf_counter() - started, str(exc), started, depth)
        raise
    else:
        report_stage(package, name, time.perf_counter() - started, None, started, depth)
    finally:
        _depth.reset(token)


class StageCollector(logging.Handler):
    """Handler that gathers every ``blade.stage`` report into one list."""

    def __init__(self) -> None:
        super().__init__(level=STAGE)
        self.reports: list[StageReport] = []

    def emit(self, record: logging.LogRecord) -> None:
        name = getattr(record, "stage_name", None)
        seconds = getattr(record, "stage_seconds", None)
        if not isinstance(name, str):
            return
        if isinstance(seconds, bool) or not isinstance(seconds, (int, float)):
            return
        error = getattr(record, "stage_error", None)
        depth = getattr(record, "stage_depth", 0)
        started = getattr(record, "stage_started", 0.0)
        package = getattr(record, "stage_package", "")
        self.reports.append(StageReport(
            str(package), name, float(seconds),
            None if error in (None, "") else str(error),
            int(depth) if isinstance(depth, int) else 0,
            float(started) if isinstance(started, (int, float)) and not isinstance(started, bool)
            else 0.0,
        ))

    def take(self) -> list[StageReport]:
        """Return the reports gathered so far and clear the list."""
        found = self.reports
        self.reports = []
        return found


@contextmanager
def collecting() -> Iterator[StageCollector]:
    """Attach a collector to ``blade.stage`` for the duration of the block."""
    log = logging.getLogger(LOGGER_NAME)
    previous = log.level
    handler = StageCollector()
    log.setLevel(STAGE)
    log.addHandler(handler)
    try:
        yield handler
    finally:
        log.removeHandler(handler)
        log.setLevel(previous)


def fold_reports(reports: Sequence[StageReport]) -> dict[str, object]:
    """One value per stage name.

    A nested report (depth above 0, or a span contained in another span) is
    dropped, so a solve inside buckling is not also the static stage. Spans of
    the same name that do not contain each other are added. An error keeps the
    seconds and the error text.
    """
    top = [item for item in reports if item.depth <= 0]
    kept = [item for item in top if not _nested_in(item, top)]
    grouped: dict[str, list[StageReport]] = {}
    for item in kept:
        grouped.setdefault(item.stage, []).append(item)
    folded: dict[str, object] = {}
    for name, items in grouped.items():
        seconds = float(sum(item.seconds for item in items))
        errors = [item.error for item in items if item.error]
        if errors:
            folded[name] = {"seconds": seconds, "error": errors[0]}
        else:
            folded[name] = seconds
    return folded


def _nested_in(item: StageReport, reports: Sequence[StageReport]) -> bool:
    """True when another report's span contains this one.

    Equal spans keep the earlier object so a package and its caller, timing the
    same call, collapse to one record.
    """
    for other in reports:
        if other is item or not _contains(other, item):
            continue
        if not _contains(item, other):
            return True
        if other.seconds > item.seconds or (other.seconds == item.seconds and id(other) < id(item)):
            return True
    return False


def _contains(outer: StageReport, inner: StageReport) -> bool:
    end_outer = outer.started + outer.seconds
    end_inner = inner.started + inner.seconds
    return inner.started >= outer.started - 1e-9 and end_inner <= end_outer + 1e-9
