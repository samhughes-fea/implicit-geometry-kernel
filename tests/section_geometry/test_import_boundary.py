import importlib
from pathlib import Path
import ast


def test_banned_homes_absent() -> None:
    root = Path(importlib.import_module("implicit_geometry_kernel.section_geometry").__file__).resolve().parent
    for banned in ("engine", "interface", "core"):
        assert not (root / banned).exists(), f"banned home present: {banned}/"


def test_production_does_not_import_material_utilities() -> None:
    root = Path(importlib.import_module("implicit_geometry_kernel.section_geometry").__file__).resolve().parent
    skip = {"tests", "visualisation", "examples", "docs", "graphify-out", "__pycache__"}
    hits: list[str] = []
    for path in root.rglob("*.py"):
        if any(part in skip for part in path.parts):
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names.extend(a.name for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                names.append(node.module)
            for name in names:
                if "material_utilities" in name or "gauge_graph" in name:
                    hits.append(f"{path.relative_to(root).as_posix()}: {name}")
    assert not hits, "section_geometry must not import material_utilities:\n" + "\n".join(hits)
