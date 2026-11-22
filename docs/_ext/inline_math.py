"""Accept ``\\( … \\)`` inline maths in the Markdown sources.

MyST reads ``\\(`` as an escaped parenthesis and drops the backslash, so pages
written with ``\\( … \\)`` reached the site as plain text. This hook rewrites
``\\( … \\)`` (on one line, outside fenced code) to ``$ … $`` before parsing, so
existing pages render while they are converted. New pages use ``$ … $``, which
GitHub renders as well.
"""

from __future__ import annotations

import re

_INLINE = re.compile(r"\\\((.+?)\\\)")
_FENCE = re.compile(r"^\s*(```+|~~~+)(.*)$")


def _rewrite(text: str) -> str:
    out, stack = [], []  # stack of "code" / "directive" fences
    for line in text.split("\n"):
        m = _FENCE.match(line)
        if m:
            info = m.group(2).strip()
            if stack and not info:
                stack.pop()
            else:
                stack.append("directive" if info.startswith("{") and not info.startswith("{code") else "code")
            out.append(line)
            continue
        out.append(line if "code" in stack else _INLINE.sub(r"$\1$", line))
    return "\n".join(out)


def _source_read(app, docname, source):
    if app.env.doc2path(docname).suffix == ".md":
        source[0] = _rewrite(source[0])


def setup(app):
    app.connect("source-read", _source_read)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
