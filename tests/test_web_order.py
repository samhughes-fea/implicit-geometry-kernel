"""IEA 15 web anchors stay in chord order through the web window (BG-037)."""

import os
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _catalogue():
    env = os.environ.get("REFERENCE_BLADES_DIR")
    roots = ([Path(env)] if env else []) + [ROOT / "reference_blades", ROOT.parent / "reference-blades"]
    for path in roots:
        if (path / "wind" / "iea15" / "IEA-15-240-RWT.yaml").is_file():
            return path
    return None


def _yaml():
    root = _catalogue()
    return root / "wind" / "iea15" / "IEA-15-240-RWT.yaml"


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_iea15_cap_stays_outside_the_webs_near_the_tip():
    import runpy

    ns = runpy.run_path(str(ROOT / "examples" / "blade" / "08_reference_blades" / "reference_blades.py"))
    _name, geo, _z = ns["load_windio"](_yaml())
    # The arc grid's last web sample is near 112.2 m. Past it the section still
    # places the arcs, and the anchors have to follow or the cap passes the web.
    for z in np.linspace(112.0, 114.0, 21):
        feet = geo._web_positions_fn(float(z))
        anchors = dict(geo.region_graph.evaluate_anchors(float(z)))
        if len(feet) < 2:
            continue
        mids = [0.5 * (a + b) for a, b in feet]
        assert anchors["web_0"] == pytest.approx(mids[0], abs=1e-6)
        assert anchors["web_1"] == pytest.approx(mids[1], abs=1e-6)
        assert anchors["web_0"] < anchors["web_1"]
        if "cap_lo" in anchors:
            assert anchors["cap_lo"] < anchors["web_0"]
        if "cap_hi" in anchors:
            assert anchors["web_1"] < anchors["cap_hi"]
    at_tip = dict(geo.region_graph.evaluate_anchors(113.3))
    assert "cap_lo" in at_tip and "cap_hi" in at_tip
    assert at_tip["cap_lo"] < at_tip["web_0"] < at_tip["web_1"] < at_tip["cap_hi"]
