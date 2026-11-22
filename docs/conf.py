"""Sphinx project for the existing Markdown tree plus a generated API reference.

Build: ``pip install -e ".[docs]"`` then ``sphinx-build -W --keep-going -b html docs docs/_build/html``
or ``make -C docs html``. The manual is the Markdown in ``docs/``, organised as four
manuals (U use, R reference theory, V validation, D development; see
``docs/conventions/DOCUMENTATION_STANDARDS.md``); ``docs/api/`` documents each
package's ``__all__`` with autodoc, so the reference is the docstrings and nothing
is copied. Citations come from ``docs/references.bib`` (sphinxcontrib-bibtex). CI builds it with warnings as errors and
``docs.yml`` publishes it to GitHub Pages from ``main``.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent / "_ext"))
_init = ast.parse((ROOT / "src" / "implicit_geometry_kernel" / "__init__.py").read_text())
release = next(n.value.value for n in _init.body if isinstance(n, ast.Assign)
               and any(getattr(t, "id", "") == "__version__" for t in n.targets))
version = release

project = "implicit-geometry-kernel"
copyright = "Sam Hughes"
author = "Sam Hughes"

extensions = [
    "myst_parser",
    "sphinxcontrib.mermaid",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxcontrib.bibtex",
    "repo_links",
    "inline_math",
]

source_suffix = {".md": "markdown"}

# One bibliography for the U, R, V and D manuals; author-year citations,
# written {cite:t}`key` (textual) or {cite:p}`key` (parenthetical).
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "plain"
bibtex_reference_style = "author_year"

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
html_title = f"implicit-geometry-kernel {release}"
html_short_title = "IGK"
html_copy_source = False
html_show_sourcelink = False
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_favicon = "_static/logo.svg"
html_theme_options = {
    "light_logo": "logo.svg",
    "dark_logo": "logo-dark.svg",
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "source_repository": "https://github.com/samhughes-fea/implicit-geometry-kernel/",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#041E42",
        "color-brand-content": "#007288",
        "font-stack": "Inter, system-ui, -apple-system, Segoe UI, Roboto, sans-serif",
    },
    "dark_css_variables": {
        "color-brand-primary": "#9DC3E6",
        "color-brand-content": "#5CC2D6",
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/samhughes-fea/implicit-geometry-kernel",
            "html": '<svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>',
            "class": "",
        },
    ],
}
copybutton_prompt_text = r">>> |\$ "
copybutton_prompt_is_regexp = True

# Relative links to examples/ and src/ are for GitHub; they are not Sphinx docs.
suppress_warnings = [
    "myst.xref_missing",
    "myst.header",
    "toc.not_included",
]
