"""
visualisation.plot
==================
Matplotlib-based visualisation helpers for section geometry.

Thin re-export surface over :mod:`plot_sdf`, :mod:`plot_section`, and
:mod:`figure_io`. Medial-axis / gradient plotters were removed from the
public API (G3).
"""

from .figure_io import (
    DEFAULT_DPI,
    PRESENTATION_DPI,
    add_output_and_presentation_args,
    apply_presentation_style,
    demo_output_dir,
    resolve_save_dpi,
    save_figure,
)
from .plot_sdf import plot_sdf_field
from .plot_section import (
    _DEFAULT_COLOR,
    _component_color,
    apply_structural_section_axes,
    fit_square_section_axes,
    plot_section,
    plot_section_from_mesh,
    plot_section_hybrid,
    plot_section_with_quadtree,
)
from .plot_topology import plot_shell_midline_strip_topology
from .plot_subcomponents import draw_section_sdf_panel, taxonomy_plot_caption

__all__ = [
    "DEFAULT_DPI",
    "PRESENTATION_DPI",
    "add_output_and_presentation_args",
    "apply_presentation_style",
    "demo_output_dir",
    "resolve_save_dpi",
    "save_figure",
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
]
