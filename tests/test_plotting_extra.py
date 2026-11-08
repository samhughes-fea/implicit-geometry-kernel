"""Canonical plots extra."""

from __future__ import annotations


def test_canonical_plot_section_import():
    from implicit_geometry_kernel.plotting.plot import plot_section

    assert callable(plot_section)
