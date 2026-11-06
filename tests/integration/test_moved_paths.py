"""The former integration.quadtree and integration.quadrature paths keep working."""

import importlib

PKG = "implicit_geometry_kernel.integration."


def test_quadtree_package_keeps_old_names():
    qt = importlib.import_module(PKG + "quadtree")
    for sub in ("tree", "refine", "refine_multi", "soa", "marching"):
        mod = importlib.import_module(PKG + "quadtree." + sub)
        names = [n for n in vars(mod) if getattr(getattr(mod, n), "__module__", None) == mod.__name__]
        for name in names:
            assert getattr(qt, name) is getattr(mod, name), name


def test_quadrature_alias():
    old = importlib.import_module(PKG + "quadrature")
    for sub in ("moments", "quadrature"):
        mod = importlib.import_module(PKG + "quadtree." + sub)
        for name in ("section_properties_from_quadtree", "label_fill_polygons", "_rectangle_moments"):
            if hasattr(mod, name) and getattr(mod, name).__module__ == mod.__name__:
                assert getattr(old, name) is getattr(mod, name)
