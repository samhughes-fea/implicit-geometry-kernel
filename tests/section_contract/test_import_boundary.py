import importlib
from pathlib import Path


def test_contract_homes_exist() -> None:
    root = Path(importlib.import_module("implicit_geometry_kernel.section_contract").__file__).resolve().parent
    assert (root / "api").is_dir()
    assert (root / "types").is_dir()
