"""Architecture gates for the G-row ownership matrix."""

from __future__ import annotations
import importlib

import ast
from pathlib import Path

import pytest

ROOT = Path(importlib.import_module("implicit_geometry_kernel.section").__file__).resolve().parent

BANNED_TOP_LEVEL = ("engine", "interface", "core")

# Skip-layer edges forbidden by G_ROW_OWNERSHIP_MATRIX.
FORBIDDEN_IMPORT_PREFIXES: dict[str, tuple[str, ...]] = {
    "layout": (
        "implicit_geometry_kernel.section.recipes",
        "implicit_geometry_kernel.section.constructive",
        "implicit_geometry_kernel.section.junctions",
        "implicit_geometry_kernel.section.assembly.components",
    ),
    "junctions": (
        "implicit_geometry_kernel.sdf",
        "implicit_geometry_kernel.section.constructive",
        "implicit_geometry_kernel.section.recipes",
        "implicit_geometry_kernel.section.layout",
    ),
    "cells": (
        "implicit_geometry_kernel.sdf",
        "implicit_geometry_kernel.section.constructive",
        "implicit_geometry_kernel.section.recipes",
    ),
    "midlines": (
        "implicit_geometry_kernel.section.constructive",
        "implicit_geometry_kernel.sdf",
        "implicit_geometry_kernel.section.props",
    ),
}


def test_banned_software_layer_homes_absent() -> None:
    for name in BANNED_TOP_LEVEL:
        assert not (ROOT / name).exists(), f"banned home still present: {name}/"


def test_api_and_recipe_homes_exist() -> None:
    assert (ROOT / "api").is_dir()
    assert (ROOT / "recipes").is_dir()
    assert (ROOT / "assembly").is_dir()
    assert (ROOT / "junctions").is_dir()
    assert (ROOT / "cells").is_dir()
    assert (ROOT / "props").is_dir()
    assert (Path(__file__).resolve().parents[2] / "docs" / "section" / "reference" / "G_ROW_OWNERSHIP_MATRIX.md").is_file()


def test_assembly_does_not_export_sdf_ir_components() -> None:
    from implicit_geometry_kernel.section import assembly

    banned = {
        "OuterSkin",
        "SparCap",
        "ShearWebComponent",
        "SandwichCore",
        "LEInsert",
        "TEInsert",
    }
    exported = set(getattr(assembly, "__all__", ()))
    assert banned.isdisjoint(exported), f"assembly still exports {banned & exported}"


def test_sections_components_dir_absent() -> None:
    assert not (ROOT / "assembly" / "components").exists()
    assert not (ROOT / "sections").exists()


_REMOVED_SHIM_PACKAGES = (
    "sections",
    "junction_preparation",
    "cell_wall_loops",
    "props_emit",
)


def test_deprecated_section_packages_absent() -> None:
    for name in _REMOVED_SHIM_PACKAGES:
        assert not (ROOT / name).exists(), f"deprecated package still present: {name}/"
    assert not (ROOT / "recipes" / "sdf_ir").exists()
    kernel = ROOT.parent
    assert not (kernel / "blade_geometry.py").exists()
    assert not (kernel / "section_geometry").exists()
    assert not (kernel / "frames" / "section_frames").exists()
    assert not (ROOT / "tagged.py").exists()
    assert not (ROOT / "eval_cache.py").exists()
    assert not (ROOT / "offset_math.py").exists()



def test_profiling_lives_outside_sdf_kernel() -> None:
    assert (ROOT.parent / "profiling.py").is_file()
    assert not (ROOT.parent / "sdf" / "profiling.py").exists()


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
