"""Former section subpackage paths are gone as of 0.6.0."""

from pathlib import Path

import implicit_geometry_kernel.section as section

ROOT = Path(section.__file__).resolve().parent
REMOVED = ("api", "placement", "constructive", "props", "recipes")


def test_moved_section_packages_absent() -> None:
    for name in REMOVED:
        assert not (ROOT / name).exists(), name
    assert not (ROOT / "midlines" / "offset_math.py").exists()
