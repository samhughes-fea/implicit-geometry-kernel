"""Run every Python block of docs/tutorial.md in order, so the tutorial stays correct."""

from __future__ import annotations

import re
from pathlib import Path

TUTORIAL = Path(__file__).resolve().parents[1] / "docs" / "tutorial.md"


def test_tutorial_blocks_run():
    blocks = re.findall(r"```python\n(.*?)```", TUTORIAL.read_text(encoding="utf-8"), re.S)
    assert len(blocks) >= 5
    namespace: dict = {}
    for block in blocks:
        exec(compile(block, str(TUTORIAL), "exec"), namespace)
    assert sorted({w.kind for w in namespace["cut"].walls}) == ["cap", "skin", "web"]
    assert set(namespace["mesh"].surface) >= {"skin", "web", "tip_closure"}
