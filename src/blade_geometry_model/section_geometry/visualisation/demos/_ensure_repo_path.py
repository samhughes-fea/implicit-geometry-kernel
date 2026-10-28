"""Ensure the repository root is on ``sys.path`` when examples are run as scripts."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def ensure_repo_path() -> Path:
    """Insert the checkout root (parent of ``blade_geometry_model/``) first on ``sys.path``.

    Raises
    ------
    RuntimeError
        If no ancestor of this file contains ``blade_geometry_model/__init__.py``.
    """
    for p in Path(__file__).resolve().parents:
        if (p / "blade_geometry_model" / "__init__.py").is_file():
            s = str(p)
            if s not in sys.path:
                sys.path.insert(0, s)
            _maybe_enable_plot_debug(p)
            return p
    raise RuntimeError(
        "Could not find the blade_geometry_model package. Clone the repository and run this "
        "script from inside it, or install the project in editable mode "
        "(from the repo root: pip install -e .)."
    )


def _maybe_enable_plot_debug(repo_root: Path) -> None:
    """Enable optional plot-debug tracing for demo scripts."""
    dbg_path = Path(__file__).resolve().parent / "_plot_debug.py"
    if not dbg_path.is_file():
        return
    try:
        spec = importlib.util.spec_from_file_location("_plot_debug", dbg_path)
        if spec is None or spec.loader is None:
            return
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        fn = getattr(mod, "activate_plot_debug_if_enabled", None)
        if callable(fn):
            fn(repo_root)
    except Exception as exc:
        print(f"[plot-debug] setup skipped: {exc}")
