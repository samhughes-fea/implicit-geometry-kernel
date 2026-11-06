"""Old sdf module paths stay importable as aliases of the new tree until 0.6.0."""

import importlib

import pytest

MOVED = {
    "field": "nodes.base",
    "nodes": "nodes.tagged",
    "span_nodes": "nodes.span",
    "lofted": "nodes.lofted",
    "primitives": "kernels.primitives",
    "_polygon_numba": "kernels.polygon_numba",
    "grid": "kernels.grid",
    "csg": "kernels.csg",
    "transforms": "kernels.transforms",
    "airfoil": "shapes.airfoil",
    "naca": "shapes.naca",
}

PKG = "implicit_geometry_kernel.sdf."


@pytest.mark.parametrize("old,new", sorted(MOVED.items()))
def test_old_path_is_new_module(old, new):
    new_mod = importlib.import_module(PKG + new)
    if old == "nodes":
        # sdf.nodes is now the package; the tagged nodes are re-exported from it.
        pkg = importlib.import_module(PKG + "nodes")
        for name in getattr(new_mod, "__all__", []):
            assert getattr(pkg, name) is getattr(new_mod, name)
        return
    assert importlib.import_module(PKG + old) is new_mod


def test_ir_reexports_shapes():
    ir = importlib.import_module(PKG + "ir")
    nodes = importlib.import_module(PKG + "nodes")
    for name in ("Circle", "Union", "Offset", "eval_expr_batch", "expr_bbox"):
        assert getattr(ir, name) is getattr(nodes, name)
