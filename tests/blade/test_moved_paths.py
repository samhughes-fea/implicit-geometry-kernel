"""``blade.cut`` still exports the shell-graph names until 0.6.0."""

import importlib
import sys

import implicit_geometry_kernel.blade.cut as cut


def test_cut_reexports_shell_graph_names():
    # ``blade.shell_graph`` is also a public function on the blade package, so
    # the submodule is taken from sys.modules (loaded by ``blade.cut``).
    mod = sys.modules.get("implicit_geometry_kernel.blade.shell_graph")
    if mod is None:
        mod = importlib.import_module("implicit_geometry_kernel.blade.shell_graph")
    assert mod is not cut.shell_graph
    for name in mod.__all__:
        assert getattr(cut, name) is getattr(mod, name), name
