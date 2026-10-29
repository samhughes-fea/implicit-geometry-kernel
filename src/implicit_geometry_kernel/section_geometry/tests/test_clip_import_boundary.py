"""Production code must import clip helpers from clip_algebra, not clips."""

from __future__ import annotations

import re
from pathlib import Path

_SECTION_GEOMETRY = Path(__file__).resolve().parents[1]
_SKIP_FILES = frozenset({"clips.py", "clip_algebra.py"})

def _iter_production_py() -> list[Path]:
    out: list[Path] = []
    for path in _SECTION_GEOMETRY.rglob("*.py"):
        rel = path.relative_to(_SECTION_GEOMETRY)
        if "tests" in rel.parts or "visualisation" in rel.parts:
            continue
        if path.name in _SKIP_FILES:
            continue
        out.append(path)
    return out

def test_no_production_imports_from_clips_module():
    failures: list[str] = []
    pattern = re.compile(
        r"from\s+([\w.]*)sections\.clips\s+import\s+"
        r"|from\s+\.\.?clips\s+import\s+"
        r"|import\s+([\w.]*)sections\.clips\b"
    )
    for path in _iter_production_py():
        text = path.read_text(encoding="utf-8")
        if "clips" not in text:
            continue
        if pattern.search(text):
            failures.append(str(path.relative_to(_SECTION_GEOMETRY)))
    assert not failures, (
        "Import clip helpers from clip_algebra, not clips.py:\n"
        + "\n".join(failures)
    )
