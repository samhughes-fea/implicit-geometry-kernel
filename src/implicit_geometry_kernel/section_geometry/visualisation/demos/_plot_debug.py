"""Matplotlib debug tracing helpers for section-geometry demo scripts."""

from __future__ import annotations

import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

_ENABLED_FLAG = "_section_geometry_plot_debug_enabled"


def _env_true(name: str) -> bool:
    v = str(os.environ.get(name, "")).strip().lower()
    return v in {"1", "true", "yes", "on", "y"}


def _json_default(obj: Any):
    if isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return {
            "shape": list(obj.shape),
            "dtype": str(obj.dtype),
        }
    return str(obj)


def _summarize_numeric_array(arr: np.ndarray) -> dict[str, Any]:
    a = np.asarray(arr)
    out: dict[str, Any] = {
        "shape": list(a.shape),
        "dtype": str(a.dtype),
        "size": int(a.size),
    }
    if a.size == 0:
        return out
    if np.issubdtype(a.dtype, np.number):
        finite = np.isfinite(a)
        out["finite_fraction"] = float(np.mean(finite))
        if np.any(finite):
            af = a[finite]
            out["min"] = float(np.min(af))
            out["max"] = float(np.max(af))
            out["mean"] = float(np.mean(af))
    return out


def _summarize_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return {"kind": "ndarray", **_summarize_numeric_array(value)}
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return {"kind": type(value).__name__, "len": 0}
        try:
            arr = np.asarray(value)
            if arr.ndim >= 1 and arr.size > 0 and np.issubdtype(arr.dtype, np.number):
                return {
                    "kind": type(value).__name__,
                    "len": len(value),
                    **_summarize_numeric_array(arr),
                }
        except Exception:
            pass
        return {"kind": type(value).__name__, "len": len(value)}
    if isinstance(value, (int, float, bool, str)):
        return value
    if value is None:
        return None
    return {"kind": type(value).__name__}


def _summarize_args(args: tuple[Any, ...], kwargs: dict[str, Any]) -> dict[str, Any]:
    out = {
        "n_args": int(len(args)),
        "args": [_summarize_value(a) for a in args[:4]],
        "kwargs": {str(k): _summarize_value(v) for k, v in list(kwargs.items())[:8]},
    }
    if len(args) > 4:
        out["args_truncated"] = True
    if len(kwargs) > 8:
        out["kwargs_truncated"] = True
    return out


def _resolve_log_path(repo_root: Path) -> Path:
    custom_dir = os.environ.get("SECTION_GEOMETRY_PLOT_DEBUG_DIR", "").strip()
    if custom_dir:
        log_dir = Path(custom_dir)
    else:
        log_dir = (
            repo_root
            / "blade_structure_utilities"
            / "geometry_utilities"
            / "section_geometry"
            / "visualisation"
            / "output"
            / "_plot_debug"
        )
    script = Path(sys.argv[0]).stem if sys.argv and sys.argv[0] else "interactive"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    pid = os.getpid()
    run_tag = os.environ.get("SECTION_GEOMETRY_PLOT_DEBUG_TAG", "").strip()
    suffix = f"_{run_tag}" if run_tag else ""
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir / f"{script}_{stamp}_pid{pid}{suffix}.jsonl"


def activate_plot_debug_if_enabled(repo_root: Path) -> Path | None:
    """Enable matplotlib plot-call tracing when `SECTION_GEOMETRY_PLOT_DEBUG=1`."""
    if not _env_true("SECTION_GEOMETRY_PLOT_DEBUG"):
        return None

    import matplotlib

    if bool(getattr(matplotlib, _ENABLED_FLAG, False)):
        return None

    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    log_path = _resolve_log_path(repo_root)

    def emit(event: dict[str, Any]) -> None:
        row = {
            "ts_utc": datetime.now(timezone.utc).isoformat(),
            "event": event,
        }
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(row, default=_json_default) + "\n")

    def wrap_axes(method_name: str) -> None:
        original = getattr(Axes, method_name, None)
        if original is None:
            return
        if bool(getattr(original, "_section_plot_debug_wrapped", False)):
            return

        def wrapped(self, *args, **kwargs):
            t0 = time.perf_counter()
            out = original(self, *args, **kwargs)
            t_ms = (time.perf_counter() - t0) * 1000.0
            try:
                emit(
                    {
                        "type": "axes_call",
                        "method": method_name,
                        "duration_ms": float(t_ms),
                        "axes_title": str(self.get_title() or ""),
                        "line_count": int(len(getattr(self, "lines", []))),
                        "collection_count": int(len(getattr(self, "collections", []))),
                        "image_count": int(len(getattr(self, "images", []))),
                        "patch_count": int(len(getattr(self, "patches", []))),
                        "payload": _summarize_args(args, kwargs),
                    }
                )
            except Exception as exc:
                emit({"type": "debug_error", "method": method_name, "error": str(exc)})
            return out

        wrapped._section_plot_debug_wrapped = True  # type: ignore[attr-defined]
        setattr(Axes, method_name, wrapped)

    methods = [
        "plot",
        "scatter",
        "bar",
        "hist",
        "boxplot",
        "contour",
        "contourf",
        "pcolormesh",
        "imshow",
        "fill",
        "fill_between",
    ]
    for name in methods:
        wrap_axes(name)

    savefig_original = Figure.savefig

    def savefig_wrapped(self, fname, *args, **kwargs):
        t0 = time.perf_counter()
        out = savefig_original(self, fname, *args, **kwargs)
        t_ms = (time.perf_counter() - t0) * 1000.0
        emit(
            {
                "type": "figure_save",
                "path": str(fname),
                "duration_ms": float(t_ms),
                "n_axes": int(len(getattr(self, "axes", []))),
                "axes": [
                    {
                        "title": str(ax.get_title() or ""),
                        "xlim": [float(ax.get_xlim()[0]), float(ax.get_xlim()[1])],
                        "ylim": [float(ax.get_ylim()[0]), float(ax.get_ylim()[1])],
                        "line_count": int(len(getattr(ax, "lines", []))),
                        "collection_count": int(len(getattr(ax, "collections", []))),
                        "image_count": int(len(getattr(ax, "images", []))),
                        "patch_count": int(len(getattr(ax, "patches", []))),
                    }
                    for ax in getattr(self, "axes", [])
                ],
            }
        )
        return out

    Figure.savefig = savefig_wrapped
    setattr(matplotlib, _ENABLED_FLAG, True)

    emit(
        {
            "type": "debug_start",
            "script": str(sys.argv[0] if sys.argv else ""),
            "python": sys.version,
            "matplotlib_backend": str(matplotlib.get_backend()),
            "log_path": str(log_path),
        }
    )
    print(f"[plot-debug] enabled -> {log_path}")
    return log_path
