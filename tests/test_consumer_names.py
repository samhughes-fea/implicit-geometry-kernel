"""The kernel names no consumer packages in src/ or tests/."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BANNED = re.compile(
    r"shell[-_]fem|anisotropic[-_]section[-_]model|composite_laminate_model|"
    r"blade[-_]structure|blade_precompute|blade_structure_utilities|"
    r"material_utilities|laminate_constitutive_model",
    re.IGNORECASE,
)
ALLOW = {
    ROOT / "tests" / "test_consumer_names.py",
    ROOT / "tests" / "test_import_boundary.py",
    ROOT / "tests" / "section" / "test_import_boundary.py",
}


def _iter_text_files() -> list[Path]:
    out: list[Path] = []
    for folder in (ROOT / "src", ROOT / "tests"):
        for path in folder.rglob("*"):
            if path.suffix.lower() not in {".py", ".md"}:
                continue
            if path.resolve() in {p.resolve() for p in ALLOW}:
                continue
            if "__pycache__" in path.parts:
                continue
            out.append(path)
    return out


def test_no_consumer_package_names_in_src_or_tests() -> None:
    hits: list[str] = []
    for path in _iter_text_files():
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if BANNED.search(line):
                rel = path.relative_to(ROOT).as_posix()
                hits.append(f"{rel}:{i}:{line.strip()}")
    assert not hits, "consumer package names still present:\n" + "\n".join(hits)
