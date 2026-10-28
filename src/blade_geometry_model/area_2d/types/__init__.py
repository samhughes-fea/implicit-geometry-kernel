"""Compatibility proxy for historical ``types.py`` imports.

This package-level surface exists for worker-contract compliance while preserving
existing ``from .types import ...`` imports in parent package modules.
"""

from __future__ import annotations

import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

_parent_pkg = __name__.rsplit(".", 1)[0]
_types_path = Path(__file__).resolve().parent.parent / "types.py"
_spec = spec_from_file_location(f"{_parent_pkg}._legacy_types", _types_path)
if _spec is None or _spec.loader is None:
    raise ImportError(f"Unable to load legacy types module at {_types_path}.")

_legacy = module_from_spec(_spec)
sys.modules[_spec.name] = _legacy
_spec.loader.exec_module(_legacy)

__all__ = ["PhiSampler", "FieldMesh2D", "AreaMesh2D"]
globals().update({name: getattr(_legacy, name) for name in __all__})
