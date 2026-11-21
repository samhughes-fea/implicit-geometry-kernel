"""Fail when a relative Markdown link in the repository points at nothing.

    python scripts/check_links.py

Checks every tracked ``*.md`` file. A link is ``[text](target)`` with a relative
target; ``http(s)``, ``mailto`` and in-page ``#anchors`` are skipped, and the
``#anchor`` part of a file link is ignored. Targets resolve against the linking
file's folder and must be a tracked file or a folder that holds one.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r"\]\(([^)#\s]+)(?:#[^)]*)?\)")
SKIP = ("http://", "https://", "mailto:")


def tracked() -> list[str]:
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True)
    return out.stdout.split()


def broken(files: list[str]) -> list[tuple[str, str]]:
    known = set(files)
    for f in files:
        d = os.path.dirname(f)
        while d:
            known.add(d)
            d = os.path.dirname(d)
    bad = []
    for f in files:
        if not f.endswith(".md"):
            continue
        text = (ROOT / f).read_text(encoding="utf-8", errors="replace")
        for m in LINK.finditer(text):
            target = m.group(1)
            if target.startswith(SKIP):
                continue
            path = os.path.normpath(os.path.join(os.path.dirname(f), target)).rstrip("/")
            if path not in known:
                bad.append((f, target))
    return bad


def main() -> int:
    bad = broken(tracked())
    for f, target in bad:
        print(f"{f}: {target}")
    print(f"{len(bad)} broken relative links")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
