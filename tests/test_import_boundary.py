"""Import rules for the kernel.

The kernel is the base of the suite:

    implicit_geometry_kernel  <-  composite_laminate_model  <-  anisotropic_section_model

so it never imports the laminate or the section model. Inside the kernel, only
``plotting`` uses Matplotlib, and ``sdf`` stands on its own (it may be split out as a
general implicit-geometry library), so it imports nothing from ``section`` or ``blade``.
"""

import ast
import importlib
import warnings
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parents[1] / "src" / "implicit_geometry_kernel"
FORBIDDEN = {"composite_laminate_model", "anisotropic_section_model", "laminate_modelling_package",
             "blade_geometry_model"}


def _imports(path: Path):
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            yield node.module
        elif isinstance(node, ast.Import):
            yield from (a.name for a in node.names)


def _library_files(exclude=()):
    for path in sorted(PACKAGE.rglob("*.py")):
        rel = path.relative_to(PACKAGE)
        if rel.parts[0] not in exclude:
            yield path, rel


def test_no_imports_of_the_models_built_on_it():
    for path, rel in _library_files():
        for module in _imports(path):
            assert module.split(".")[0] not in FORBIDDEN, (rel, module)


def test_only_plotting_uses_matplotlib():
    for path, rel in _library_files(exclude=("plotting",)):
        for module in _imports(path):
            if rel.as_posix() == "sdf/grid.py" and module.startswith("matplotlib"):
                continue  # contour fallback inside a function when contourpy is missing
            assert not module.startswith("matplotlib"), (rel, module)
            assert not module.startswith("implicit_geometry_kernel.plotting"), (rel, module)


def test_sdf_stands_alone():
    for path, rel in _library_files():
        if rel.parts[0] != "sdf":
            continue
        for module in _imports(path):
            if module.startswith("implicit_geometry_kernel."):
                top = module.split(".")[1]
                assert top in {"sdf", "profiling"}, (rel, module)


@pytest.mark.parametrize("module, name", [
    ("implicit_geometry_kernel.blade.airfoil", "Station"),
    ("implicit_geometry_kernel.section.layout.system_layout", "SYSTEM_TYPE_KEYS"),
])
def test_deprecated_aliases_still_import(module, name):
    import sys

    sys.modules.pop(module, None)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mod = importlib.import_module(module)
    assert hasattr(mod, name)
    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
