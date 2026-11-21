"""Forwarding for the deprecated flat ``mesh.*`` module paths (0.9.2 layout)."""

from __future__ import annotations

import importlib
import warnings
from typing import Any, Callable


def forward(old: str, new: str) -> Callable[[str], Any]:
    """Module ``__getattr__`` that warns once per read and returns ``new``'s attribute."""

    def __getattr__(name: str) -> Any:
        if name.startswith("__"):
            raise AttributeError(name)
        target = importlib.import_module(new)
        try:
            value = getattr(target, name)
        except AttributeError:
            raise AttributeError(f"module {old!r} has no attribute {name!r}") from None
        warnings.warn(
            f"{old} is deprecated; import {name} from {new} (or implicit_geometry_kernel.mesh). "
            "The flat path is removed in 0.9.6.",
            DeprecationWarning,
            stacklevel=2,
        )
        return value

    return __getattr__
