"""The geometry kernel is the base of the split:

    implicit_geometry_kernel  <-  composite_laminate_model  <-  anisotropic_section_model

so it never imports the laminate or the section model.
"""

import ast
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1] / "src" / "implicit_geometry_kernel"
FORBIDDEN = {"composite_laminate_model", "anisotropic_section_model", "laminate_modelling_package",
             "blade_geometry_model"}


def test_no_imports_of_the_models_built_on_it():
    for path in sorted(PACKAGE.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                modules = [node.module]
            elif isinstance(node, ast.Import):
                modules = [a.name for a in node.names]
            else:
                continue
            for module in modules:
                assert module.split(".")[0] not in FORBIDDEN, (path.relative_to(PACKAGE), module)
