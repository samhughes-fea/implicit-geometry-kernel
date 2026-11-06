"""Locked section-plane structural frame identifiers (S / S′ / B).

Source of truth for readable frame names. Letters remain short prose/math
shortcuts after the long name is introduced once.
"""

from __future__ import annotations

# S — zero-twist station; origin on the pitch axis
ZERO_TWIST_PITCH_AXIS_STATION_FRAME = "zero_twist_pitch_axis_station_frame"

# S′ — zero-twist station; origin on the centroid axis
ZERO_TWIST_CENTROID_AXIS_STATION_FRAME = "zero_twist_centroid_axis_station_frame"

# B — twisted blade-installed attitude; origin on the pitch axis
TWISTED_PITCH_AXIS_BLADE_FRAME = "twisted_pitch_axis_blade_frame"

FRAME_ID_BY_LETTER: dict[str, str] = {
    "S": ZERO_TWIST_PITCH_AXIS_STATION_FRAME,
    "S'": ZERO_TWIST_CENTROID_AXIS_STATION_FRAME,
    "S′": ZERO_TWIST_CENTROID_AXIS_STATION_FRAME,
    "B": TWISTED_PITCH_AXIS_BLADE_FRAME,
}

LETTER_BY_FRAME_ID: dict[str, str] = {
    ZERO_TWIST_PITCH_AXIS_STATION_FRAME: "S",
    ZERO_TWIST_CENTROID_AXIS_STATION_FRAME: "S′",
    TWISTED_PITCH_AXIS_BLADE_FRAME: "B",
}

__all__ = [
    "ZERO_TWIST_PITCH_AXIS_STATION_FRAME",
    "ZERO_TWIST_CENTROID_AXIS_STATION_FRAME",
    "TWISTED_PITCH_AXIS_BLADE_FRAME",
    "FRAME_ID_BY_LETTER",
    "LETTER_BY_FRAME_ID",
]
