"""The former frames.tensors module path keeps every name it had."""

import importlib


def test_tensors_package_keeps_old_names():
    pkg = importlib.import_module("implicit_geometry_kernel.frames.tensors")
    for sub in ("inertia", "moments", "thickness", "forces"):
        mod = importlib.import_module(f"implicit_geometry_kernel.frames.tensors.{sub}")
        for name, obj in vars(mod).items():
            if getattr(obj, "__module__", None) == mod.__name__:
                assert getattr(pkg, name) is obj, name
