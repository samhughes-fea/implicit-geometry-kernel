"""Ensure the repository root is on ``sys.path`` when examples are run as scripts."""

from __future__ import annotations

import sys
from pathlib import Path


def ensure_repo_path() -> Path:
    """Insert the checkout root (parent of ``implicit_geometry_kernel/``) first on ``sys.path``.

    Raises
    ------
    RuntimeError
        If no ancestor of this file contains ``implicit_geometry_kernel/__init__.py``.
    """
    for p in Path(__file__).resolve().parents:
        if (p / "implicit_geometry_kernel" / "__init__.py").is_file():
            s = str(p)
            if s not in sys.path:
                sys.path.insert(0, s)
            return p
    raise RuntimeError(
        "Could not find the implicit_geometry_kernel package. Clone the repository and run this "
        "script from inside it, or install the project in editable mode "
        "(from the repo root: pip install -e .)."
    )
