"""Former section subpackage paths stay importable as aliases until 0.6.0."""

import importlib
import pkgutil

import pytest

PKG = "implicit_geometry_kernel.section."
MOVED = {"props": "properties", "recipes": "walls"}


def _cases():
    for old, new in MOVED.items():
        pkg = importlib.import_module(PKG + new)
        for info in pkgutil.iter_modules(pkg.__path__):
            yield old, new, info.name


@pytest.mark.parametrize("old,new,mod", list(_cases()))
def test_old_submodule_is_new_module(old, new, mod):
    assert importlib.import_module(f"{PKG}{old}.{mod}") is importlib.import_module(f"{PKG}{new}.{mod}")


@pytest.mark.parametrize("old,new", sorted(MOVED.items()))
def test_old_package_reexports(old, new):
    o, n = importlib.import_module(PKG + old), importlib.import_module(PKG + new)
    names = getattr(n, "__all__", [k for k in vars(n) if not k.startswith("_")])
    for name in names:
        assert getattr(o, name) is getattr(n, name), name


CONSTRUCTIVE = {
    "emitters": "boundary.emitters",
    "labels": "boundary.labels",
    "properties": "properties.polygons",
}


@pytest.mark.parametrize("old,new", sorted(CONSTRUCTIVE.items()))
def test_constructive_modules_alias(old, new):
    assert importlib.import_module(f"{PKG}constructive.{old}") is importlib.import_module(PKG + new)


PLACEMENT = {
    "anchoring": "layout.anchoring",
    "collision": "layout.collision",
    "handoff": "layout.handoff",
}


@pytest.mark.parametrize("old,new", sorted(PLACEMENT.items()))
def test_placement_modules_alias(old, new):
    assert importlib.import_module(f"{PKG}placement.{old}") is importlib.import_module(PKG + new)


def test_placement_package_reexports_layout_modules():
    placement = importlib.import_module(PKG + "placement")
    anchoring = importlib.import_module(PKG + "layout.anchoring")
    assert placement.AnchoredWallGeometry is anchoring.AnchoredWallGeometry
    assert placement.anchored_wall_geometry is anchoring.anchored_wall_geometry


def test_constructive_export_reexports_walls_and_cores():
    export = importlib.import_module(PKG + "constructive.export")
    walls = importlib.import_module(PKG + "boundary.walls")
    cores = importlib.import_module(PKG + "boundary.cores")
    assert export.emit_wall_solid_polygons is walls.emit_wall_solid_polygons
    assert export.emit_core_polygons is cores.emit_core_polygons
    assert export.emit_constructive_polygons is cores.emit_constructive_polygons
