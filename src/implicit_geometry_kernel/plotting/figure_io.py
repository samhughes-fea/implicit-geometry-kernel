"""Shared figure I/O for section_geometry visualisation demos.

Demos write their figures under ``examples/demos/figures/<demo_stem>/`` (passed in as
``root``) via :func:`demo_output_dir` and :func:`save_figure`. The figures are committed
and embedded in the demo READMEs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib as mpl

DEFAULT_DPI = 150
PRESENTATION_DPI = 220


def demo_output_dir(
    demo_stem: str,
    *,
    repo_visualisation_root: Path | None = None,
) -> Path:
    """Return ``<root>/figures/<demo_stem>/`` (created on save, not here).

    Demos pass ``examples/demos`` as the root; without one the figures land in the
    current working directory, never inside the installed package.
    """
    root = Path(repo_visualisation_root) if repo_visualisation_root is not None else Path.cwd()
    return (root / "figures" / str(demo_stem)).resolve()


def _chrome_rc(*, font_size: int, title_size: int, label_size: int, tick_size: int,
               legend_size: int, figure_title_size: int) -> dict[str, Any]:
    from .palette import INK, MUTED

    return {
        "font.size": font_size,
        "axes.titlesize": title_size,
        "axes.labelsize": label_size,
        "xtick.labelsize": tick_size,
        "ytick.labelsize": tick_size,
        "legend.fontsize": legend_size,
        "figure.titlesize": figure_title_size,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.8,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "grid.color": "#d7dce3",
        "grid.linewidth": 0.6,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }


def apply_doc_style() -> None:
    """Quiet chrome and deck ink for README / docs figures (not slide scale)."""
    mpl.rcParams.update(
        _chrome_rc(
            font_size=10,
            title_size=11,
            label_size=10,
            tick_size=9,
            legend_size=9,
            figure_title_size=12,
        )
    )


def apply_presentation_style(presentation: bool) -> None:
    """Apply slide chrome when ``presentation`` is true, otherwise doc style.

    Fonts, but also the chrome: a figure headed for a slide wants a lighter
    frame and a quieter grid than one headed for a docs page, and it wants the
    deck's ink colour rather than pure black.
    """
    if not presentation:
        apply_doc_style()
        return
    mpl.rcParams.update(
        _chrome_rc(
            font_size=12,
            title_size=14,
            label_size=12,
            tick_size=11,
            legend_size=10,
            figure_title_size=16,
        )
    )


def save_figure(
    fig: Any,
    path: Path | str,
    *,
    dpi: int | None = None,
    presentation: bool = False,
    **savefig_kwargs: Any,
) -> Path:
    """Save ``fig`` with white face and tight bbox; create parent dirs."""
    if presentation:
        apply_presentation_style(True)
    else:
        apply_doc_style()
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if dpi is None:
        dpi = PRESENTATION_DPI if presentation else DEFAULT_DPI
    kwargs = {
        "dpi": int(dpi),
        "facecolor": "white",
        "bbox_inches": "tight",
    }
    kwargs.update(savefig_kwargs)
    fig.savefig(out, **kwargs)
    return out


def add_output_and_presentation_args(
    parser: argparse.ArgumentParser,
    *,
    default_stem: str,
    root: Path | None = None,
) -> None:
    """Add ``--out-dir`` and ``--presentation``; the default is ``<root>/figures/<default_stem>/``."""
    default_out = demo_output_dir(default_stem, repo_visualisation_root=root)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=default_out,
        help=f"Output directory (default: figures/{default_stem}/).",
    )
    parser.add_argument(
        "--presentation",
        action="store_true",
        help="Larger fonts and higher DPI for slides/decks.",
    )


def resolve_save_dpi(*, presentation: bool, dpi: int | None = None) -> int:
    """Pick DPI from explicit override, presentation flag, or default."""
    if dpi is not None:
        return int(dpi)
    return PRESENTATION_DPI if presentation else DEFAULT_DPI


__all__ = [
    "DEFAULT_DPI",
    "PRESENTATION_DPI",
    "add_output_and_presentation_args",
    "apply_doc_style",
    "apply_presentation_style",
    "demo_output_dir",
    "resolve_save_dpi",
    "save_figure",
]
