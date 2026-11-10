"""Sphinx project for the existing Markdown tree. Not a second manual.

Build: ``pip install -e ".[docs]"`` then ``sphinx-build -b html docs docs/_build/html``
or ``make -C docs html``. Do not enable autodoc: the glossary is the name index.
GitHub Pages / Read the Docs are not wired here.
"""

from __future__ import annotations

project = "implicit-geometry-kernel"
copyright = "Sam Hughes"
author = "Sam Hughes"

extensions = [
    "myst_parser",
    "sphinxcontrib.mermaid",
]

# Markdown only. No sphinx.ext.autodoc / autosummary.
source_suffix = {".md": "markdown"}
root_doc = "index"
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**/*.aux",
    "**/*.log",
    "**/*.toc",
]

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
]
myst_fence_as_directive = ["mermaid"]
myst_heading_anchors = 3

html_theme = "alabaster"
html_title = "implicit-geometry-kernel"
html_copy_source = True

# Relative links to examples/ and src/ are for GitHub; they are not Sphinx docs.
suppress_warnings = [
    "myst.xref_missing",
    "myst.header",
    "toc.not_included",
]
