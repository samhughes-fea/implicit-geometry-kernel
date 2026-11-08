"""Structural family labels on a catalogue token (letter X)."""

from __future__ import annotations

from typing import Literal

StructuralFamily = Literal["A", "B", "C", "D"]


def parse_structural_family(value: str | StructuralFamily) -> str:
    s = str(value).strip().upper()
    if s not in ("A", "B", "C", "D"):
        raise ValueError(
            f"structural_family must be one of A|B|C|D, got {value!r}"
        )
    return s


__all__ = ["StructuralFamily", "parse_structural_family"]
