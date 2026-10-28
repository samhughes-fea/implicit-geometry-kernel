"""Extend G-row architecture gates for the final hardening pass."""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _iter_prod_py() -> list[Path]:
    skip = {"tests", "visualisation", "examples", "docs", "graphify-out", "__pycache__"}
    out: list[Path] = []
    for p in ROOT.rglob("*.py"):
        if any(part in skip for part in p.parts):
            continue
        out.append(p)
    return out


def test_sections_components_dir_absent() -> None:
    assert not (ROOT / "sections" / "components").exists()


def test_no_from_component_outside_tests() -> None:
    violations: list[str] = []
    for path in _iter_prod_py():
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name == "from_component":
                violations.append(path.relative_to(ROOT).as_posix())
            if isinstance(node, ast.AsyncFunctionDef) and node.name == "from_component":
                violations.append(path.relative_to(ROOT).as_posix())
    assert not violations, "from_component still defined in:\n" + "\n".join(violations)


def test_no_cap_from_components_symbol() -> None:
    hits: list[str] = []
    for path in _iter_prod_py():
        text = path.read_text(encoding="utf-8")
        if "_cap_from_components" in text:
            hits.append(path.relative_to(ROOT).as_posix())
    assert not hits, "_cap_from_components still present in:\n" + "\n".join(hits)


def test_no_legacy_vocab_in_public_all() -> None:
    banned = {
        "emit_laminate_polygons",
        "is_laminate_label",
        "build_void_cell_wall_loops",
    }
    hits: list[str] = []
    for path in _iter_prod_py():
        if path.name != "__init__.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id == "__all__":
                        if isinstance(node.value, (ast.List, ast.Tuple)):
                            for elt in node.value.elts:
                                if isinstance(elt, ast.Constant) and elt.value in banned:
                                    hits.append(
                                        f"{path.relative_to(ROOT).as_posix()}: {elt.value}"
                                    )
    assert not hits, "legacy vocab still in __all__:\n" + "\n".join(hits)


def test_no_laminate_identifiers_in_production() -> None:
    """Wall solids own geometric thickness; laminate vocab belongs to material workers.

    Private one-release aliases live only in ``layout/thickness_limits.py``.
    """
    import re

    allow = { (ROOT / "layout" / "thickness_limits.py").resolve() }
    pat = re.compile(r"\blaminate\b", re.IGNORECASE)
    hits: list[str] = []
    for path in _iter_prod_py():
        if path.resolve() in allow:
            continue
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if pat.search(line):
                hits.append(f"{path.relative_to(ROOT).as_posix()}:{i}:{line.strip()}")
    assert not hits, "laminate still in production section_geometry:\n" + "\n".join(hits)


_G_ROW_PUBLIC_ALL = frozenset(
    {
        "SectionGeometryResult",
        "SectionGeometryRunner",
        "SectionPropertiesMeta",
        "StationLayoutResolution",
        "SectionLayout",
        "SectionMidlineHost",
        "build_section_geometry",
        "build_section_with_thicknesses",
        "section_geometry_prep_key",
        "AirfoilSDF",
        "MultiCellSection",
        "resolve_station_layout",
        "resolve_system_type",
        "SYSTEM_TYPE_KEYS",
        "SystemLayoutSpec",
        "layout_without_spar_caps",
        "is_core_label",
        "MIN_REALISTIC_SKIN_WALL_THICKNESS_M",
        "MIN_REALISTIC_SPAR_WALL_THICKNESS_M",
        "MIN_REALISTIC_WEB_WALL_THICKNESS_M",
        "clamp_skin_thickness_m",
        "clamp_spar_thickness_m",
        "clamp_web_thickness_m",
    }
)

_PROPS_DEBUG_BARREL_BAN = frozenset(
    {
        "props_backend_from_env",
        "quadtree_spec_for_airfoil",
        "build_shared_section_quadtree",
    }
)


def _package_all() -> set[str]:
    init = ROOT / "__init__.py"
    tree = ast.parse(init.read_text(encoding="utf-8"), filename=str(init))
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "__all__":
                    if isinstance(node.value, (ast.List, ast.Tuple)):
                        return {
                            elt.value
                            for elt in node.value.elts
                            if isinstance(elt, ast.Constant) and isinstance(elt.value, str)
                        }
    raise AssertionError("package __all__ not found")


def test_g_row_public_barrel_allowlist() -> None:
    actual = _package_all()
    assert actual == _G_ROW_PUBLIC_ALL, (
        "package __all__ drifted from G-row / L4-EVAL surface:\n"
        f"  missing: {sorted(_G_ROW_PUBLIC_ALL - actual)}\n"
        f"  extra: {sorted(actual - _G_ROW_PUBLIC_ALL)}"
    )


def test_props_debug_not_in_public_all() -> None:
    actual = _package_all()
    hits = sorted(actual & _PROPS_DEBUG_BARREL_BAN)
    assert not hits, f"G-PROPS debug symbols still in package __all__: {hits}"


def test_no_cap_recipes_from_assembled_symbol() -> None:
    hits: list[str] = []
    for path in _iter_prod_py():
        if "_cap_recipes_from_assembled" in path.read_text(encoding="utf-8"):
            hits.append(path.relative_to(ROOT).as_posix())
    assert not hits, "_cap_recipes_from_assembled still present in:\n" + "\n".join(hits)


def test_sections_barrel_does_not_eager_import_multicell() -> None:
    init = ROOT / "sections" / "__init__.py"
    tree = ast.parse(init.read_text(encoding="utf-8"), filename=str(init))
    for node in tree.body:
        if isinstance(node, (ast.ImportFrom,)):
            if node.level == 1 and node.module == "multicell":
                raise AssertionError(
                    "sections/__init__.py eagerly imports multicell "
                    "(use lazy __getattr__ to keep the assemble DAG acyclic)"
                )


def test_sections_relative_import_dag_acyclic() -> None:
    """Relative imports among sections/*.py must form a DAG."""
    from collections import defaultdict, deque

    sec = ROOT / "sections"
    graph: dict[str, set[str]] = defaultdict(set)
    modules = {p.stem for p in sec.glob("*.py") if p.name != "__init__.py"}
    for path in sec.glob("*.py"):
        if path.name == "__init__.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        src = path.stem
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom) or node.level != 1:
                continue
            if node.module is None:
                # ``from . import clips`` → edge to clips, not package __init__
                for alias in node.names:
                    if alias.name in modules:
                        graph[src].add(alias.name)
            elif node.module.split(".", 1)[0] in modules:
                graph[src].add(node.module.split(".", 1)[0])

    indeg = {m: 0 for m in modules}
    for src, dsts in graph.items():
        for d in dsts:
            if d in indeg:
                indeg[d] += 1
            else:
                indeg[d] = 1
                modules.add(d)
    q = deque([m for m, n in indeg.items() if n == 0])
    seen = 0
    while q:
        u = q.popleft()
        seen += 1
        for v in graph.get(u, ()):
            indeg[v] -= 1
            if indeg[v] == 0:
                q.append(v)
    assert seen == len(indeg), (
        "cyclic relative imports in sections/: "
        + ", ".join(sorted(m for m, n in indeg.items() if n > 0))
    )
