"""Compare local chord fractions vs root-anchored web positions (drops at tip).

Uses layout ``5D-F`` on a short tip chord so outboard anchored webs drop.

Run from the repository root::

    python src/blade_geometry_model/section_geometry/examples/04_spanwise_web_anchors.py
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
if not _erp.is_file():
    raise RuntimeError(f"Missing helper {_erp.name}; run from the repository checkout.")
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    SectionGeometryRunner,
    resolve_system_type,
)

ROOT_CHORD_M = 0.2409
TIP_CHORD_M = 0.0878
T_SKIN = 0.001


def _report(title: str, resolution) -> None:
    print(f"\n=== {title} ===")
    print(f"nominal={resolution.nominal_layout.key} effective={resolution.layout.key}")
    print(f"survivors_m={tuple(round(x, 4) for x in resolution.web_x_positions_m)}")
    if resolution.dropped:
        for d in resolution.dropped:
            print(f"  dropped x_m={d.x_m:.4f} reason={d.reason}")
    else:
        print("  dropped=(none)")
    print(f"stepped_down={resolution.stepped_down}")


def main() -> None:
    layout = resolve_system_type("5D-F")
    tip = AirfoilSDF.from_naca("0012", n_points=200, chord=TIP_CHORD_M)
    anchored = tuple(float(f) * ROOT_CHORD_M for f in layout.web_chord_fracs)

    legacy = (
        SectionGeometryRunner(
            layout=layout,
            airfoil=tip,
            skin_thickness=T_SKIN,
            spar_cap_thickness=0.0,
            shear_web_thickness=0.0,
            # omit web_x_positions_m → frac × local_chord
        )
        .resolve()
        .get_layout_resolution()
    )
    _report("legacy local fractions (no drops for taper)", legacy)

    anchored_res = (
        SectionGeometryRunner(
            layout=layout,
            airfoil=tip,
            skin_thickness=T_SKIN,
            spar_cap_thickness=0.0,
            shear_web_thickness=0.0,
            web_x_positions_m=anchored,
        )
        .resolve()
        .get_layout_resolution()
    )
    _report("root-anchored positions (drops outboard webs)", anchored_res)


if __name__ == "__main__":
    main()
