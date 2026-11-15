"""Kernel stages report on the shared logger without changing the profile sink."""

from __future__ import annotations

import pytest

from implicit_geometry_kernel.profiling import (
    LOGGER_NAME,
    StageCollector,
    StageReport,
    collecting,
    fold_reports,
    stage,
    timed_stage,
)


def test_a_stage_records_seconds_and_an_error_keeps_its_time() -> None:
    with collecting() as sink:
        with stage("implicit_geometry_kernel", "kept"):
            pass
        with pytest.raises(RuntimeError):
            with stage("implicit_geometry_kernel", "failed"):
                raise RuntimeError("stage failed")
    folded = fold_reports(sink.take())
    assert folded["kept"] >= 0.0
    failed = folded["failed"]
    assert isinstance(failed, dict)
    assert failed["seconds"] >= 0.0
    assert failed["error"] == "stage failed"


def test_timed_stage_reports_and_still_fills_its_sink() -> None:
    rows: dict[str, float] = {}
    with collecting() as sink:
        with timed_stage("shell-prep", sink=rows):
            pass
    assert rows["shell-prep"] >= 0.0
    folded = fold_reports(sink.take())
    assert folded["shell-prep"] >= 0.0


def test_a_nested_span_keeps_measures_and_is_dropped_from_the_total() -> None:
    with collecting() as sink:
        with stage("implicit_geometry_kernel", "shell_mesh") as outer:
            outer.set("elements", 12)
            outer.set("edge_length", 0.05)
            with stage("implicit_geometry_kernel", "inner"):
                pass
    reports = sink.take()
    assert reports[0].depth == 1
    assert reports[1].depth == 0
    assert reports[1].measures["elements"] == 12
    assert reports[1].measures["edge_length"] == 0.05
    assert reports[1].cpu_seconds >= 0.0
    assert set(fold_reports(reports)) == {"shell_mesh"}


def test_public_contract() -> None:
    import dataclasses
    import logging

    assert LOGGER_NAME == "blade.stage"
    assert [item.name for item in dataclasses.fields(StageReport)] == [
        "package", "stage", "seconds", "error", "depth", "started",
        "cpu_seconds", "measures"]
    assert callable(stage) and callable(fold_reports)
    assert issubclass(StageCollector, logging.Handler)
