"""The integration.quadtree package keeps its split-module names."""

from pathlib import Path

import importlib

import implicit_geometry_kernel.integration as integration

PKG = "implicit_geometry_kernel.integration."


def test_quadtree_package_keeps_old_names():
    qt = importlib.import_module(PKG + "quadtree")
    for sub in ("tree", "soa", "marching"):
        mod = importlib.import_module(PKG + "quadtree." + sub)
        names = [n for n in vars(mod) if getattr(getattr(mod, n), "__module__", None) == mod.__name__]
        for name in names:
            assert getattr(qt, name) is getattr(mod, name), name


def test_quadrature_alias_removed():
    root = Path(integration.__file__).resolve().parent
    assert not (root / "quadrature.py").exists()
