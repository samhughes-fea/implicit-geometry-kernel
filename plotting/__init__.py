"""Deprecated import path; use ``implicit_geometry_kernel.plotting``. Removed in 0.6.0."""

from __future__ import annotations

import sys
import warnings

warnings.warn(
    "import plotting is deprecated; use implicit_geometry_kernel.plotting. Removed in 0.6.0.",
    DeprecationWarning,
    stacklevel=1,
)

import implicit_geometry_kernel.plotting as _moved

_prefix = "implicit_geometry_kernel.plotting"
for _name, _mod in list(sys.modules.items()):
    if _name == _prefix or _name.startswith(_prefix + "."):
        sys.modules["plotting" + _name[len(_prefix) :]] = _mod
sys.modules[__name__] = _moved
