"""Golden numbers and public API stay put unless a PR updates them on purpose.

Regenerate with ``python scripts/snapshots.py --update`` and say why in CHANGELOG.md.
"""

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("snapshots", ROOT / "scripts" / "snapshots.py")
snapshots = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(snapshots)

SECTIONS = json.loads(snapshots.SECTIONS_FILE.read_text())


@pytest.mark.parametrize("case", sorted(SECTIONS))
def test_section_golden(case):
    key, twist = case.split("@")
    moved = snapshots.diff(SECTIONS[case], snapshots.section_case(key, float(twist)))
    assert not moved, "\n".join(moved[:20])


def test_blade_golden():
    moved = snapshots.diff(json.loads(snapshots.BLADE_FILE.read_text()), snapshots.compute_blade())
    assert not moved, "\n".join(moved[:20])


def test_public_api_unchanged():
    moved = snapshots.diff(json.loads(snapshots.API_FILE.read_text()), snapshots.compute_api())
    removed_or_changed = [m for m in moved if not m.endswith(": new")]
    assert not removed_or_changed, (
        "public API changed; deprecate first (see CHANGELOG.md), then update the snapshot:\n"
        + "\n".join(removed_or_changed[:20])
    )


def test_every_layout_is_covered():
    from implicit_geometry_kernel import SYSTEM_TYPE_LAYOUT_KEYS

    covered = {case.split("@")[0] for case in SECTIONS}
    assert covered == set(SYSTEM_TYPE_LAYOUT_KEYS)
