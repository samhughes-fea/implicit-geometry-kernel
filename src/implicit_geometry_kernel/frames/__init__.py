"""
Geometric coordinate frames R / B / S / S′ and the shell basis C, with point and
tensor transforms. The laminate frames E / M keep their tags here; their transforms
live in laminate-constitutive-model.

Ownership home: this package (R/B/S/S′); frame ids and tags cover the shell basis C.
Short tags: ``frame_tags`` (``FrameTag`` / ``FRAME_*`` / ``require_frame_tag``).
Demos: ``examples/demos/frames/``.

Prefer this barrel for public imports. Deep home imports are allowed for
maintainers; do not reimplement twist / NACA pivot / flap↔My/Mz outside these
homes.
"""

from implicit_geometry_kernel.frames.frame_tags import FRAME_B, FRAME_C, FRAME_E, FRAME_M, FRAME_R, FRAME_S, FRAME_SPRIME, FRAME_TAGS, FrameTag, require_frame_tag
from implicit_geometry_kernel.frames.facade import CoordinateFrameTransforms
from implicit_geometry_kernel.frames.frame_ids import FRAME_ID_BY_LETTER, LETTER_BY_FRAME_ID, TWISTED_PITCH_AXIS_BLADE_FRAME, ZERO_TWIST_CENTROID_AXIS_STATION_FRAME, ZERO_TWIST_PITCH_AXIS_STATION_FRAME
from implicit_geometry_kernel.frames.geometry import PIVOT_FRACTION_DEFAULT, S_to_naca, blade_to_chord, chord_to_blade, naca_to_S, naca_to_twisted_pitch_axis_blade_frame, naca_to_zero_twist_pitch_axis_station_frame, rotate_B_to_S, rotate_S_to_B, rotate_S_to_B_swapped_order, rotate_twisted_pitch_axis_blade_frame_to_zero_twist_pitch_axis_station_frame, rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame, rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame_swapped_order, twisted_pitch_axis_blade_frame_to_naca, zero_twist_pitch_axis_station_frame_to_naca
from implicit_geometry_kernel.frames.tensors import structural_product_moment
from implicit_geometry_kernel.frames.point_transforms import rotate_points_2d, transform_points_S_to_B, transform_points_B_to_S, transform_points_S_to_Sprime, transform_points_Sprime_to_S, transform_points_aerofoil_to_S, rR_from_zL_scalars

__all__ = [
    "FRAME_B",
    "FRAME_C",
    "FRAME_E",
    "FRAME_ID_BY_LETTER",
    "FRAME_M",
    "FRAME_R",
    "FRAME_S",
    "FRAME_SPRIME",
    "FRAME_TAGS",
    "FrameTag",
    "LETTER_BY_FRAME_ID",
    "PIVOT_FRACTION_DEFAULT",
    "TWISTED_PITCH_AXIS_BLADE_FRAME",
    "ZERO_TWIST_CENTROID_AXIS_STATION_FRAME",
    "ZERO_TWIST_PITCH_AXIS_STATION_FRAME",
    "CoordinateFrameTransforms",
    "S_to_naca",
    "blade_to_chord",
    "chord_to_blade",
    "naca_to_S",
    "naca_to_twisted_pitch_axis_blade_frame",
    "naca_to_zero_twist_pitch_axis_station_frame",
    "require_frame_tag",
    "rotate_B_to_S",
    "rotate_S_to_B",
    "rotate_S_to_B_swapped_order",
    "rotate_twisted_pitch_axis_blade_frame_to_zero_twist_pitch_axis_station_frame",
    "rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame",
    "rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame_swapped_order",
    "twisted_pitch_axis_blade_frame_to_naca",
    "zero_twist_pitch_axis_station_frame_to_naca",
    "structural_product_moment",
    "rotate_points_2d",
    "transform_points_S_to_B",
    "transform_points_B_to_S",
    "transform_points_S_to_Sprime",
    "transform_points_Sprime_to_S",
    "transform_points_aerofoil_to_S",
    "rR_from_zL_scalars",
]
