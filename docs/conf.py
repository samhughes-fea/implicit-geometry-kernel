"""Sphinx project for the existing Markdown tree plus a generated API reference.

Build: ``pip install -e ".[docs]"`` then ``sphinx-build -W --keep-going -b html docs docs/_build/html``
or ``make -C docs html``. The manual is the Markdown already in the repository;
``docs/api/`` documents each package's ``__all__`` with autodoc, so the reference
is the docstrings and nothing is copied. CI builds it with warnings as errors and
``docs.yml`` publishes it to GitHub Pages from ``main``.
"""

from __future__ import annotations

project = "implicit-geometry-kernel"
copyright = "Sam Hughes"
author = "Sam Hughes"

extensions = [
    "myst_parser",
    "sphinxcontrib.mermaid",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

source_suffix = {".md": "markdown"}

# API reference: each package's ``__all__``, in source order, signatures without
# repeated type hints (they are in the signature already).
autodoc_default_options = {"members": True, "member-order": "bysource", "show-inheritance": True}
autodoc_typehints = "signature"
autodoc_preserve_defaults = True
autodoc_mock_imports = ["gmsh", "cadquery", "contourpy"]
napoleon_google_docstring = False
napoleon_numpy_docstring = True
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

html_theme = "furo"
html_title = "implicit-geometry-kernel"
html_copy_source = True

# Relative links to examples/ and src/ are for GitHub; they are not Sphinx docs.
suppress_warnings = [
    "myst.xref_missing",
    "myst.header",
    "toc.not_included",
]
