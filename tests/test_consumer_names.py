"""The kernel names no consumer packages in src/ or tests/."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BANNED = re.compile(
    r"shell[-_]fem|anisotropic[-_]section[-_]model|composite_laminate_model|"
    r"blade[-_]structure|blade_precompute|blade_structure_utilities|"
    r"material_utilities|laminate_constitutive_model|tidalshell",
    re.IGNORECASE,
)
ALLOW = {
    ROOT / "tests" / "test_consumer_names.py",
    ROOT / "tests" / "test_import_boundary.py",
    ROOT / "tests" / "section" / "test_import_boundary.py",
}


def _iter_text_files() -> list[Path]:
    out: list[Path] = []
    for folder in (ROOT / "src", ROOT / "tests", ROOT / "examples"):
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


# Private consumers are never named anywhere in this public repository, docs and
# changelog included. Public sibling packages may be named in docs and history.
PRIVATE = re.compile(r"tidalshell", re.IGNORECASE)
TEXT = {".py", ".md", ".toml", ".yml", ".yaml", ".txt", ".cfg", ".rst"}


def test_no_private_consumer_names_anywhere() -> None:
    hits: list[str] = []
    for path in sorted(ROOT.rglob("*")):
        rel = path.relative_to(ROOT)
        if path.suffix.lower() not in TEXT or not path.is_file():
            continue
        if rel.parts[0] in {".git", "build", "dist", "reference_blades"} or "__pycache__" in rel.parts:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        for i, line in enumerate(path.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
            if PRIVATE.search(line):
                hits.append(f"{rel.as_posix()}:{i}:{line.strip()}")
    assert not hits, "private consumer names present:\n" + "\n".join(hits)
