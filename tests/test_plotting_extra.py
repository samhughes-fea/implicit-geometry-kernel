"""Canonical plots extra and root alias until 0.6.0."""

from __future__ import annotations

import sys
import warnings


def test_canonical_plot_section_import():
    from implicit_geometry_kernel.plotting.plot import plot_section

    assert callable(plot_section)


def test_root_plotting_alias_warns_until_060():
    sys.modules.pop("plotting", None)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        import plotting  # noqa: F401
        from plotting.plot import plot_section as aliased

    assert any(
        issubclass(item.category, DeprecationWarning) and "0.6.0" in str(item.message) for item in caught
    )
    from implicit_geometry_kernel.plotting.plot import plot_section

    assert aliased is plot_section
    assert sys.modules["plotting"] is sys.modules["implicit_geometry_kernel.plotting"]
