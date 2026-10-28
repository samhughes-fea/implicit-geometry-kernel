"""Section frames B / S / S′ (+ R span grids).

See ``README.md`` for mathematics and Exact I/O. Prefer the package barrel
``coordinate_frames`` for public imports.
"""

from .facade import CoordinateFrameTransforms
from .frame_ids import (
    FRAME_ID_BY_LETTER,
    LETTER_BY_FRAME_ID,
    TWISTED_PITCH_AXIS_BLADE_FRAME,
    ZERO_TWIST_CENTROID_AXIS_STATION_FRAME,
    ZERO_TWIST_PITCH_AXIS_STATION_FRAME,
)
from .geometry import (
    PIVOT_FRACTION_DEFAULT,
    S_to_naca,
    blade_to_chord,
    chord_to_blade,
    naca_to_S,
    naca_to_twisted_pitch_axis_blade_frame,
    naca_to_zero_twist_pitch_axis_station_frame,
    rotate_B_to_S,
    rotate_S_to_B,
    rotate_S_to_B_swapped_order,
    rotate_twisted_pitch_axis_blade_frame_to_zero_twist_pitch_axis_station_frame,
    rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame,
    rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame_swapped_order,
    twisted_pitch_axis_blade_frame_to_naca,
    zero_twist_pitch_axis_station_frame_to_naca,
)
from .plot import (
    SECTION_PLOT_XLABEL,
    SECTION_PLOT_YLABEL,
    fit_square_section_axes,
    plot_section_polyline,
    section_normal_2d,
    section_plot_axes,
    setup_section_axes,
)
from .span import SpanFrameTransforms
from .tensors import structural_product_moment
from .point_transforms import (
    rotate_points_2d,
    transform_points_S_to_B,
    transform_points_B_to_S,
    transform_points_S_to_Sprime,
    transform_points_Sprime_to_S,
    transform_points_aerofoil_to_S,
    rR_from_zL_scalars,
)

__all__ = [
    "CoordinateFrameTransforms",
    "SpanFrameTransforms",
    "structural_product_moment",
    "rotate_points_2d",
    "transform_points_S_to_B",
    "transform_points_B_to_S",
    "transform_points_S_to_Sprime",
    "transform_points_Sprime_to_S",
    "transform_points_aerofoil_to_S",
    "rR_from_zL_scalars",
    "FRAME_ID_BY_LETTER",
    "LETTER_BY_FRAME_ID",
    "TWISTED_PITCH_AXIS_BLADE_FRAME",
    "ZERO_TWIST_CENTROID_AXIS_STATION_FRAME",
    "ZERO_TWIST_PITCH_AXIS_STATION_FRAME",
    "PIVOT_FRACTION_DEFAULT",
    "S_to_naca",
    "blade_to_chord",
    "chord_to_blade",
    "naca_to_S",
    "naca_to_twisted_pitch_axis_blade_frame",
    "naca_to_zero_twist_pitch_axis_station_frame",
    "rotate_B_to_S",
    "rotate_S_to_B",
    "rotate_S_to_B_swapped_order",
    "rotate_twisted_pitch_axis_blade_frame_to_zero_twist_pitch_axis_station_frame",
    "rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame",
    "rotate_zero_twist_pitch_axis_station_frame_to_twisted_pitch_axis_blade_frame_swapped_order",
    "twisted_pitch_axis_blade_frame_to_naca",
    "zero_twist_pitch_axis_station_frame_to_naca",
    "SECTION_PLOT_XLABEL",
    "SECTION_PLOT_YLABEL",
    "fit_square_section_axes",
    "plot_section_polyline",
    "section_normal_2d",
    "section_plot_axes",
    "setup_section_axes",
]
