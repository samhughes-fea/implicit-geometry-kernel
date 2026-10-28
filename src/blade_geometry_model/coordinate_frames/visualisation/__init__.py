"""Operational visualisation suite for coordinate_frames (demos + output).

Plot helpers and demos; frame mathematics live in ``section_frames`` /
``shell_frames``. Re-exports the section plot contract for demo convenience.
"""

from blade_geometry_model.coordinate_frames.section_frames.plot import (
    SECTION_PLOT_XLABEL,
    SECTION_PLOT_YLABEL,
    fit_square_section_axes,
    plot_section_polyline,
    section_plot_axes,
    setup_section_axes,
)

__all__ = [
    "SECTION_PLOT_XLABEL",
    "SECTION_PLOT_YLABEL",
    "fit_square_section_axes",
    "plot_section_polyline",
    "section_plot_axes",
    "setup_section_axes",
]
