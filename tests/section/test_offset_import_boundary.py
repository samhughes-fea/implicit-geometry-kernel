"""No Shapely in the repository. GEOS lives in ``implicit_geometry_kernel.geos``.

Shoelace sums numpy rings. Plotting uses ``oriented_rings``.
"""

from __future__ import annotations

import ast
import importlib
from pathlib import Path

ROOT = Path(importlib.import_module("implicit_geometry_kernel").__file__).resolve().parents[2]
PACKAGE = Path(importlib.import_module("implicit_geometry_kernel").__file__).resolve().parent
REMOVED_ALIAS_FILES = {
    "section/constructive/helpers.py",
    "section/constructive/offset_locus.py",
    "section/midlines/offset_math.py",
    "section/boundary/geos.py",
}


def _imports(path: Path):
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            yield node.module
        elif isinstance(node, ast.Import):
            yield from (a.name for a in node.names)


def _py_files(*roots: Path):
    for root in roots:
        for path in sorted(root.rglob("*.py")):
            yield path


def test_repository_does_not_import_shapely():
    roots = (PACKAGE, ROOT / "tests", ROOT / "examples")
    bad = [
        (str(path.relative_to(ROOT)), m)
        for path in _py_files(*roots)
        for m in _imports(path)
        if m.split(".")[0] == "shapely"
    ]
    assert not bad, f"Shapely imports remain: {bad}"


def test_offset_alias_modules_removed():
    for rel in REMOVED_ALIAS_FILES:
        assert not (PACKAGE / rel).exists(), rel
