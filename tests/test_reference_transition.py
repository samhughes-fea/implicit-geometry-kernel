"""Published blades: cylinder–airfoil stations should cut after plateau collapse."""

from pathlib import Path

import pytest

from implicit_geometry_kernel.contract import InfeasibleSectionGeometry

ROOT = Path(__file__).resolve().parents[1]


def _catalogue():
    for path in (ROOT / "reference_blades", ROOT.parent / "reference-blades"):
        if (path / "nrel5mw" / "nrel5mw.yaml").is_file():
            return path
    return None


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_nrel_transition_section_at():
    import runpy

    ns = runpy.run_path(str(ROOT / "examples" / "09_reference_blades" / "reference_blades.py"))
    name, geo, _ = ns["load_windio"](_catalogue() / "nrel5mw" / "nrel5mw.yaml")
    cut = geo.section_at(3.07)
    assert len(cut.walls) >= 3
    kinds = {w.kind for w in cut.walls}
    assert "skin" in kinds
