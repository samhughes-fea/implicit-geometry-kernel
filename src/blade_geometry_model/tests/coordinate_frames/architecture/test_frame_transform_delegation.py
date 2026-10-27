"""ARCH-MUST-016 — transform delegation and FrameTag field hygiene.

Reporting-mode allowlist for legitimate local spaces / dated debt. Physics
regression tests remain hard-fail. Promote allowlist shrinkage in follow-ons.
"""

from __future__ import annotations

from pathlib import Path

import pytest

_SRC = Path(__file__).resolve().parents[3]
_COORD = _SRC / "coordinate_frames"

# Packages scanned for banned hand-built transforms (the forked geometry only;
# blade-structure's mesh, load, FE and material packages did not come across).
_SCAN_ROOTS = (_SRC,)

_SKIP_DIR_NAMES = frozenset(
    {
        "tests",
        "benchmarks",
        "output",
        "__pycache__",
        "archive",
        "visualisation",
        "examples",
        "demos",
    }
)

# Canon homes may define the matrices.
_PLY_TRANSFORM_ALLOW = frozenset(
    {
        (_COORD / "shell_frames" / "ply.py").resolve(),
    }
)

# Dated debt / documented local spaces (FRAME inventory). Empty after FRAME-005 removal.
_PLY_TRANSFORM_DATED = frozenset()


def _iter_package_sources(roots: tuple[Path, ...]):
    for root in roots:
        if not root.is_dir():
            continue
        for path in root.rglob("*.py"):
            if any(part in _SKIP_DIR_NAMES for part in path.parts):
                continue
            yield path


@pytest.mark.frame_comparison
def test_no_hand_built_ply_transform_matrices_outside_canon():
    """``c * c, s * s`` 3×3 blocks are Tσ/Tε — owned by coordinate_frames."""
    offenders: list[str] = []
    allowed = _PLY_TRANSFORM_ALLOW | _PLY_TRANSFORM_DATED
    for path in _iter_package_sources(_SCAN_ROOTS):
        if path.resolve() in allowed:
            continue
        src = path.read_text(encoding="utf-8")
        if "c * c, s * s" in src or "c*c, s*s" in src:
            offenders.append(str(path.relative_to(_SRC)))
    assert offenders == [], (
        "ply stress/strain transforms must come from coordinate_frames "
        f"(ARCH-MUST-016 / FRAME-011): {offenders}"
    )


@pytest.mark.frame_comparison
def test_legacy_frame_apis_removed():
    """FRAME-004 / FRAME-005: dated wrappers are gone from live surfaces."""
    midline = _SRC / "section_geometry" / "midlines" / "shell_midline_export.py"
    midline_src = midline.read_text(encoding="utf-8")
    assert "def rotate_chord_to_blade" not in midline_src
    assert '"rotate_chord_to_blade"' not in midline_src
    assert "c * p[:, 0] - s * p[:, 1]" not in midline_src


@pytest.mark.frame_comparison
def test_frames_inventory_exists():
    inv = _COORD / "FRAMES_INVENTORY.md"
    assert inv.is_file(), "FRAMES_INVENTORY.md missing (FRAME-010)"
    text = inv.read_text(encoding="utf-8")
    assert "ARCH-MUST-016" in text
    assert "FRAME-001" in text
