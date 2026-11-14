"""Kernel stages report on the shared logger without changing the profile sink."""

from __future__ import annotations

import pytest

from implicit_geometry_kernel.profiling import timed_stage
from implicit_geometry_kernel.stage_log import collecting, fold_reports, stage


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
