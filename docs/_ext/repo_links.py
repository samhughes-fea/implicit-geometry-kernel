"""Turn Markdown links to repository files outside ``docs/`` into GitHub links.

The manual is written for GitHub, so pages link to ``../examples/...`` or
``../src/...``. Sphinx cannot resolve those as documents; this extension points
them at the file on GitHub instead of leaving a dead link.
"""

from __future__ import annotations

import posixpath
from pathlib import Path

from docutils import nodes

REPO = "https://github.com/samhughes-fea/implicit-geometry-kernel"
BRANCH = "main"


def _resolve(app, env, node, contnode):
    if node.get("refdomain") not in (None, "", "std") or node.get("reftype") not in ("myst", "doc", "ref"):
        return None
    target = str(node.get("reftarget", ""))
    if not target or "://" in target or target.startswith("#"):
        return None
    path, _, anchor = target.partition("#")
    docdir = posixpath.dirname(str(env.doc2path(node["refdoc"], False)).replace("\\", "/"))
    rel = posixpath.normpath(posixpath.join("docs", docdir, path))
    root = Path(app.confdir).parent
    if rel.startswith("..") or not (root / rel).exists():
        return None
    kind = "tree" if (root / rel).is_dir() else "blob"
    url = f"{REPO}/{kind}/{BRANCH}/{rel}" + (f"#{anchor}" if anchor else "")
    ref = nodes.reference("", "", internal=False, refuri=url)
    ref += contnode
    return ref


def setup(app):
    app.connect("missing-reference", _resolve)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
