"""Shapely has one gateway: ``section/boundary/polygons.py``.

Inside the kernel only that module and ``integration/shoelace.py`` (exact polygon
sums, a lower layer) import Shapely. Production code imports polygon and offset
helpers from ``boundary``, never from the old alias modules.
"""

from __future__ import annotations

import ast
import importlib
from pathlib import Path

PACKAGE = Path(importlib.import_module("implicit_geometry_kernel").__file__).resolve().parent
SHAPELY_HOMES = {"section/boundary/polygons.py", "integration/shoelace.py"}
ALIASES = {
    "implicit_geometry_kernel.section.constructive.helpers",
    "implicit_geometry_kernel.section.constructive.offset_locus",
    "implicit_geometry_kernel.section.midlines.offset_math",
}
ALIAS_FILES = {"section/constructive/helpers.py", "section/constructive/offset_locus.py",
               "section/midlines/offset_math.py"}


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
        if rel not in SHAPELY_HOMES and not rel.startswith("plotting/")
        for m in _imports(path)
        if m.split(".")[0] == "shapely"
    ]
    assert not bad, f"import Shapely through section.boundary: {bad}"


def test_production_code_does_not_use_the_alias_modules():
    bad = [(rel, m) for path, rel in _files() if rel not in ALIAS_FILES
           for m in _imports(path) if m in ALIASES]
    assert not bad, f"import from section.boundary instead: {bad}"


def test_aliases_point_at_the_gateway():
    boundary = importlib.import_module("implicit_geometry_kernel.section.boundary")
    for name in ALIASES:
        mod = importlib.import_module(name)
        for attr in mod.__all__:
            assert getattr(mod, attr) is getattr(boundary, attr), (name, attr)
