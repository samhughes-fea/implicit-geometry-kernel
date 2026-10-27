"""Production code must import offset APIs from offset_locus, not helpers."""

from __future__ import annotations

import re
from pathlib import Path

_SECTION_GEOMETRY = Path(__file__).resolve().parents[1]
_FORBIDDEN = frozenset(
    {"offset_polygon_inward", "airfoil_offset_ring", "airfoil_polygon"}
)
_SKIP_FILES = frozenset({"helpers.py", "offset_locus.py"})


def _iter_production_py() -> list[Path]:
    out: list[Path] = []
    for path in _SECTION_GEOMETRY.rglob("*.py"):
        rel = path.relative_to(_SECTION_GEOMETRY)
        if "tests" in rel.parts:
            continue
        if path.name in _SKIP_FILES:
            continue
        out.append(path)
    return out


def test_no_production_offset_imports_from_helpers():
    failures: list[str] = []
    pattern = re.compile(
        r"from\s+([\w.]*)constructive\.helpers\s+import\s+([^\n]+)"
        r"|from\s+\.helpers\s+import\s+([^\n]+)"
    )
    for path in _iter_production_py():
        text = path.read_text(encoding="utf-8")
        if "helpers" not in text:
            continue
        if not any(name in text for name in _FORBIDDEN):
            continue
        for m in pattern.finditer(text):
            imported = m.group(2) or m.group(3) or ""
            names = {
                part.strip().split(" as ")[0].strip() for part in imported.split(",")
            }
            bad = sorted(names & _FORBIDDEN)
            if bad:
                rel = path.relative_to(_SECTION_GEOMETRY)
                failures.append(f"{rel}: imports {bad} from helpers")
    assert not failures, (
        "Offset APIs must come from offset_locus, not helpers:\n" + "\n".join(failures)
    )
