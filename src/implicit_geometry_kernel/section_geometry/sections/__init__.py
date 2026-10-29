"""Section builders — assemble host only (G-ASSEMBLE).

SDF IR classes live under ``recipes.sdf_ir`` (private projection internals).
Public assemble symbols: ``MultiCellSection``, ``multicell_from_layout``.

``MultiCellSection`` is lazy-exported so ``clip_algebra`` / clip helpers can
load without forming an import cycle through this barrel.
"""

from __future__ import annotations

from typing import Any

from .from_layout import (
    as_section_layout,
    fixed_cap_anchor_from_layout,
    multicell_from_layout,
    multicell_from_section_layout,
)
from .labels import ComponentLabel

__all__ = [
    "ComponentLabel",
    "MultiCellSection",
    "as_section_layout",
    "fixed_cap_anchor_from_layout",
    "multicell_from_layout",
    "multicell_from_section_layout",
]


def __getattr__(name: str) -> Any:
    if name == "MultiCellSection":
        from .multicell import MultiCellSection

        return MultiCellSection
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
