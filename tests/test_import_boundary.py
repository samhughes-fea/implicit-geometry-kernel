"""Import rules for the kernel.

The kernel is the base of the suite:

    implicit_geometry_kernel  <-  shell_fem  (with laminate_constitutive_model)

and the section model reads its cuts, so it never imports the laminate, shell or section
packages. Matplotlib lives in ``implicit_geometry_kernel.plotting`` (extra ``[plots]``),
not in geometry modules.
``sdf`` stands on its own (it may be split out as a general implicit-geometry
library), so it imports nothing from ``section`` or ``blade``.
"""

import ast
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parents[1] / "src" / "implicit_geometry_kernel"
FORBIDDEN = {"laminate_constitutive_model", "shell_fem", "composite_laminate_model", "anisotropic_section_model",
             "laminate_modelling_package", "blade_geometry_model"}


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


def test_kernel_does_not_import_matplotlib_or_plotting():
    for path, rel in _library_files(exclude=("plotting",)):
        for module in _imports(path):
            assert not module.startswith("matplotlib"), (rel, module)
            assert module.split(".")[0] != "plotting", (rel, module)
            if module.startswith("implicit_geometry_kernel."):
                assert module.split(".")[1] != "plotting", (rel, module)


def test_geometry_does_not_import_mesh_or_cad_extra():
    """``mesh/`` and ``blade.cad`` are extras. ``blade/oml.py`` may lazy-import cad."""
    for path, rel in _library_files(exclude=("mesh", "plotting")):
        posix = rel.as_posix()
        for module in _imports(path):
            if module.startswith("implicit_geometry_kernel.mesh") or module == "gmsh":
                raise AssertionError((rel, module))
            if posix in {"blade/oml.py", "blade/cad.py"}:
                continue
            if module in {"implicit_geometry_kernel.blade.cad", "cadquery"}:
                raise AssertionError((rel, module))


def test_sdf_stands_alone():
    for path, rel in _library_files():
        if rel.parts[0] != "sdf":
            continue
        for module in _imports(path):
            if module.startswith("implicit_geometry_kernel."):
                top = module.split(".")[1]
                assert top in {"sdf", "profiling"}, (rel, module)


@pytest.mark.parametrize("rel", [
    "blade_geometry.py",
    "section_geometry",
    "frames/section_frames",
    "section/layout/system_layout.py",
    "blade/concept.py",
    "section/spec/factory.py",
])
def test_deprecated_top_level_aliases_removed(rel):
    assert not (PACKAGE / rel).exists(), f"0.3.0 removed path still present: {rel}"


def test_structure_sits_between_sdf_and_section():
    """``structure`` builds on ``sdf`` only; ``sdf`` never imports ``structure``."""
    for path, rel in _library_files():
        for module in _imports(path):
            if not module.startswith("implicit_geometry_kernel."):
                continue
            top = module.split(".")[1]
            if rel.parts[0] == "structure":
                assert top in {"sdf", "structure", "contract", "profiling"}, (rel, module)
            if rel.parts[0] == "sdf":
                assert top != "structure", (rel, module)


# Down-only layers: contract <- sdf, integration, frames <- structure <- catalogue <- section <- blade <- extract
_LAYER_ALLOWED = {
    "contract": {"contract", "profiling"},
    "sdf": {"sdf", "profiling"},
    "integration": {"contract", "sdf", "integration", "profiling"},
    "frames": {"contract", "sdf", "frames", "profiling"},
    "structure": {"contract", "sdf", "structure", "profiling"},
    "catalogue": {"contract", "catalogue", "profiling"},
    "section": {
        "contract", "sdf", "integration", "frames", "structure", "catalogue", "section",
        "profiling",
    },
    "blade": {
        "contract", "sdf", "integration", "frames", "structure", "catalogue", "section", "blade",
        "extract", "profiling",
    },
    "extract": {
        "contract", "sdf", "integration", "frames", "structure", "catalogue", "section", "blade",
        "extract", "profiling",
    },
    "mesh": {"contract", "blade", "extract", "mesh", "profiling"},
}


def test_package_layers_import_down_only():
    for path, rel in _library_files(exclude=("plotting",)):
        layer = rel.parts[0]
        allowed = _LAYER_ALLOWED.get(layer)
        if allowed is None:
            continue
        for module in _imports(path):
            if not module.startswith("implicit_geometry_kernel."):
                continue
            top = module.split(".")[1]
            assert top in allowed, (rel, module, allowed)
