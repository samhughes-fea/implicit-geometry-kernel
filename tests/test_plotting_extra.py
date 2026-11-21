"""Canonical plots extra."""

from __future__ import annotations


def test_canonical_plot_section_import():
    from implicit_geometry_kernel.plotting.plot import plot_section

    assert callable(plot_section)


def test_panel_colour_is_stable_and_in_the_retro_palette():
    from implicit_geometry_kernel.plotting import RETRO_PANELS, panel_colour

    roles = ("te_su", "te_pr", "cap_su", "shear_web_1", "rib:cell_0")
    colours = [panel_colour(role) for role in roles]
    assert colours == [panel_colour(role) for role in roles]
    assert all(colour in RETRO_PANELS for colour in colours)
    assert len(set(colours)) > 1
