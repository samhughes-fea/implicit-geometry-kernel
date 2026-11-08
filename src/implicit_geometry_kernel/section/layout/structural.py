"""Structural family labels and fixed-cap anchors for station assembly."""

from __future__ import annotations

from typing import Literal

from implicit_geometry_kernel.catalogue.structural import (
    StructuralFamily,
    parse_structural_family,
)

FixedCapAnchor = Literal["pitching", "max_thickness"]


def parse_fixed_cap_anchor(value: str | FixedCapAnchor) -> str:
    v = str(value).strip().lower().replace("-", "_")
    if v == "max thickness":
        v = "max_thickness"
    if v not in ("pitching", "max_thickness"):
        raise ValueError(
            f"fixed_cap_anchor must be 'pitching' or 'max_thickness', got {value!r}"
        )
    return v


__all__ = [
    "FixedCapAnchor",
    "StructuralFamily",
    "parse_fixed_cap_anchor",
    "parse_structural_family",
]
