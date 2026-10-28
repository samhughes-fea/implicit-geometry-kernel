"""Canonical short frame tags for structural containers across ``src``.

Source of truth for the R / B / S / S′ / C / E / M letter vocabulary. Long-form
section IDs live in :mod:`.section_frames.frame_ids`; transforms live in
``section_frames`` / ``shell_frames``. Containers should type ``frame`` /
``*_frame`` fields as :data:`FrameTag` and validate with
:func:`require_frame_tag` at construction or finalize boundaries.

Notes
-----
Canonical short tags are single letters (plus ``S'``). Only exact
:data:`FRAME_TAGS` members are accepted (ASCII apostrophe on S′; no alias
rewriting).

Frame roles
-----------
R
    Rotor / global reference (span metadata; not a section XYZ).
B
    Twisted pitch-axis blade frame (installed outline = S rotated by twist).
S
    Twist-agnostic station; origin on the pitch / spanwise axis.
S'
    Twist-agnostic station; origin on the centroid axis.
E
    Laminate / ply / ABD material frame.
M
    Ply-material constitutive frame.
C
    Midline / shell-panel contour frame.
"""

from __future__ import annotations

from typing import Literal

# Canonical short tags used on structural frame fields
FrameTag = Literal["R", "B", "S", "S'", "E", "M", "C"]

# Ordered catalogue of valid FrameTag values [-]
FRAME_TAGS: tuple[FrameTag, ...] = ("R", "B", "S", "S'", "E", "M", "C")

FRAME_R: FrameTag = "R"
FRAME_B: FrameTag = "B"
FRAME_S: FrameTag = "S"
FRAME_SPRIME: FrameTag = "S'"
FRAME_E: FrameTag = "E"
FRAME_M: FrameTag = "M"
FRAME_C: FrameTag = "C"


def require_frame_tag(tag: str | None) -> FrameTag:
    """Require an exact canonical :data:`FrameTag` (no alias rewriting).

    Parameters
    ----------
    tag : str
        Exact tag from :data:`FRAME_TAGS` (e.g. ``"S'"``, ``"B"``).

    Returns
    -------
    FrameTag
        The same value, typed as :data:`FrameTag`.

    Raises
    ------
    ValueError
        If ``tag`` is ``None`` or not an exact member of :data:`FRAME_TAGS`.
        Message starts with ``ERROR: frame does not exist``.
    """
    if tag is None or tag not in FRAME_TAGS:
        raise ValueError(f"ERROR: frame does not exist: {tag!r}")
    return tag  # type: ignore[return-value]


__all__ = [
    "FRAME_B",
    "FRAME_C",
    "FRAME_E",
    "FRAME_M",
    "FRAME_R",
    "FRAME_S",
    "FRAME_SPRIME",
    "FRAME_TAGS",
    "FrameTag",
    "require_frame_tag",
]
