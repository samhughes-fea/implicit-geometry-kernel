from pathlib import Path
import ast


def test_contract_homes_exist() -> None:
    root = Path(__file__).resolve().parents[1]
    assert (root / "api").is_dir()
    assert (root / "types").is_dir()
    for banned in ("engine", "interface", "core"):
        assert not (root / banned).exists(), f"banned home present: {banned}/"


def test_production_does_not_import_material_utilities() -> None:
    root = Path(__file__).resolve().parents[1]
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
