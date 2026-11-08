"""Shapely has one gateway: ``section/boundary/polygons.py``.

Shoelace sums numpy rings and must not import Shapely. Plotting uses
``oriented_rings``. Production code imports polygon helpers from ``boundary``.
"""

from __future__ import annotations

import ast
import importlib
from pathlib import Path

PACKAGE = Path(importlib.import_module("implicit_geometry_kernel").__file__).resolve().parent
SHAPELY_HOMES = {"section/boundary/polygons.py"}
REMOVED_ALIAS_FILES = {
    "section/constructive/helpers.py",
    "section/constructive/offset_locus.py",
    "section/midlines/offset_math.py",
}


def _imports(path: Path):
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            yield node.module
        elif isinstance(node, ast.Import):
            yield from (a.name for a in node.names)


def _files():
    for path in sorted(PACKAGE.rglob("*.py")):
        yield path, path.relative_to(PACKAGE).as_posix()


def test_only_the_gateway_imports_shapely():
    bad = [
        (rel, m)
        for path, rel in _files()
        if rel not in SHAPELY_HOMES
        for m in _imports(path)
        if m.split(".")[0] == "shapely"
    ]
    assert not bad, f"import Shapely through section.boundary: {bad}"


def test_offset_alias_modules_removed():
    for rel in REMOVED_ALIAS_FILES:
        assert not (PACKAGE / rel).exists(), rel
