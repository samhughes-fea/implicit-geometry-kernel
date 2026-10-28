"""Matplotlib plot helpers and visualisation demos for section_geometry."""

from .plot import (
    _DEFAULT_COLOR,
    _component_color,
    apply_structural_section_axes,
    draw_section_sdf_panel,
    fit_square_section_axes,
    plot_sdf_field,
    plot_section,
    plot_section_from_mesh,
    plot_section_hybrid,
    plot_section_with_quadtree,
    plot_shell_midline_strip_topology,
    taxonomy_plot_caption,
)
from .plot_cell_wall_loops import (
    LOOP_STYLES,
    draw_cell_wall_loops_dashed,
    draw_cell_wall_loops_filled,
    draw_component_midlines,
    plot_section_with_cell_wall_loops,
)

__all__ = [
    "plot_sdf_field",
    "apply_structural_section_axes",
    "fit_square_section_axes",
    "plot_section",
    "plot_section_from_mesh",
    "plot_section_hybrid",
    "plot_section_with_quadtree",
    "plot_shell_midline_strip_topology",
    "draw_section_sdf_panel",
    "taxonomy_plot_caption",
    "_DEFAULT_COLOR",
    "_component_color",
    "LOOP_STYLES",
    "draw_component_midlines",
    "draw_cell_wall_loops_dashed",
    "draw_cell_wall_loops_filled",
    "plot_section_with_cell_wall_loops",
]
