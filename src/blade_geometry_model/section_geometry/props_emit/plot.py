"""Backward-compatible re-export of plot helpers (now in ``visualisation``)."""

from blade_geometry_model.section_geometry.visualisation.plot import (  # noqa: F401
    _DEFAULT_COLOR,
    _component_color,
    apply_structural_section_axes,
    fit_square_section_axes,
    plot_sdf_field,
    plot_section,
    plot_section_from_mesh,
    plot_section_hybrid,
    plot_section_with_quadtree,
    plot_shell_midline_strip_topology,
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
    "_DEFAULT_COLOR",
    "_component_color",
]
