"""Architecture gates for the G-row ownership matrix."""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

BANNED_TOP_LEVEL = ("engine", "interface", "core")

# Skip-layer edges forbidden by G_ROW_OWNERSHIP_MATRIX.
FORBIDDEN_IMPORT_PREFIXES: dict[str, tuple[str, ...]] = {
    "layout": (
        "blade_structure_utilities.geometry_utilities.section_geometry.recipes",
        "blade_structure_utilities.geometry_utilities.section_geometry.constructive",
        "blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation",
        "blade_structure_utilities.geometry_utilities.section_geometry.topology",
        "blade_structure_utilities.geometry_utilities.section_geometry.sections.components",
    ),
    "junction_preparation": (
        "blade_structure_utilities.geometry_utilities.section_geometry.sdf",
        "blade_structure_utilities.geometry_utilities.section_geometry.constructive",
        "blade_structure_utilities.geometry_utilities.section_geometry.recipes",
        "blade_structure_utilities.geometry_utilities.section_geometry.layout",
    ),
    "cell_wall_loops": (
        "blade_structure_utilities.geometry_utilities.section_geometry.sdf",
        "blade_structure_utilities.geometry_utilities.section_geometry.constructive",
        "blade_structure_utilities.geometry_utilities.section_geometry.recipes",
    ),
    "midlines": (
        "blade_structure_utilities.geometry_utilities.section_geometry.constructive",
        "blade_structure_utilities.geometry_utilities.section_geometry.sdf",
        "blade_structure_utilities.geometry_utilities.section_geometry.props_emit",
    ),
}


def test_banned_software_layer_homes_absent() -> None:
    for name in BANNED_TOP_LEVEL:
        assert not (ROOT / name).exists(), f"banned home still present: {name}/"


def test_api_and_types_homes_exist() -> None:
    assert (ROOT / "api").is_dir()
    assert (ROOT / "types").is_dir()
    assert (ROOT / "recipes" / "sdf_ir").is_dir()
    assert (ROOT / "docs" / "reference" / "G_ROW_OWNERSHIP_MATRIX.md").is_file()


def test_sections_does_not_export_sdf_ir_components() -> None:
    from blade_structure_utilities.geometry_utilities.section_geometry import sections

    banned = {
        "OuterSkin",
        "SparCap",
        "ShearWeb",
        "SandwichCore",
        "LEInsert",
        "TEInsert",
    }
    exported = set(getattr(sections, "__all__", ()))
    assert banned.isdisjoint(exported), f"sections still exports {banned & exported}"


def test_sections_components_dir_absent() -> None:
    assert not (ROOT / "sections" / "components").exists()



def test_profiling_lives_outside_sdf_kernel() -> None:
    assert (ROOT / "profiling.py").is_file()


def _iter_py_files(home: str) -> list[Path]:
    root = ROOT / home
    if not root.is_dir():
        return []
    return [p for p in root.rglob("*.py") if p.is_file()]


def _module_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    mods: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                mods.add(alias.name)
        elif isinstance(node, ast.ImportFrom) and node.module:
            mods.add(node.module)
    return mods


@pytest.mark.parametrize("home,forbidden", list(FORBIDDEN_IMPORT_PREFIXES.items()))
def test_g_row_forbidden_imports(home: str, forbidden: tuple[str, ...]) -> None:
    violations: list[str] = []
    for path in _iter_py_files(home):
        for mod in _module_imports(path):
            if any(mod == p or mod.startswith(p + ".") for p in forbidden):
                rel = path.relative_to(ROOT).as_posix()
                violations.append(f"{rel} imports {mod}")
    assert not violations, "G-row skip-layer imports:\n" + "\n".join(violations)
