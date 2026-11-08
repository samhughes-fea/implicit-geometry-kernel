"""Former sdf module paths are gone as of 0.6.0."""

import ast
import importlib
from pathlib import Path

import pytest

import implicit_geometry_kernel.sdf as sdf

ROOT = Path(sdf.__file__).resolve().parent
PACKAGE = ROOT.parent
REMOVED = (
    "field.py",
    "ir.py",
    "span_nodes.py",
    "lofted.py",
    "primitives.py",
    "_polygon_numba.py",
    "grid.py",
    "csg.py",
    "transforms.py",
    "airfoil.py",
    "naca.py",
    "compile.py",
)
RETIRED_MODULES = tuple(
    f"implicit_geometry_kernel.sdf.{Path(name).stem}" for name in REMOVED
)


def test_old_sdf_modules_absent() -> None:
    for name in REMOVED:
        assert not (ROOT / name).exists(), name


@pytest.mark.parametrize("mod", RETIRED_MODULES)
def test_retired_sdf_modules_are_not_importable(mod: str) -> None:
    with pytest.raises(ImportError):
        importlib.import_module(mod)


def test_live_code_does_not_import_retired_sdf_modules() -> None:
    banned = set(RETIRED_MODULES)
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module in banned:
                rel = path.relative_to(PACKAGE)
                raise AssertionError(f"{rel} still imports {node.module}")
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name in banned:
                        rel = path.relative_to(PACKAGE)
                        raise AssertionError(f"{rel} still imports {alias.name}")
