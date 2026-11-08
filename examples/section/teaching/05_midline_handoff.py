"""Print contract midlines: the mesh hand-off from this kernel.

A 3-D shell mesh is built downstream from ``shell_graph`` polylines. Cap
footprints are ``LaminateOffsetRegion``s (constitutive OFFSET), not a second
surface and not a kinked stack midline. This script does not call a mesher.

Run from the repository root::

    python examples/section/teaching/05_midline_handoff.py
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

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.contract import LaminateOffsetRegion, ShellMidlineStrip
from implicit_geometry_kernel import shell_strips_from_graph


def main() -> None:
    stations = [
        Station(0.0, 1.0, Airfoil.naca4("2412"), pitch_axis=0.3),
        Station(1.0, 1.0, Airfoil.naca4("2412"), pitch_axis=0.3),
    ]
    cut = BladeGeometry(
        stations,
        "2D-F-LC-0R",
        0.004,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005, pivot_fraction=0.3).section_at(0.5)
    graph = cut.shell_graph()
    regions = graph.offset_regions()
    print(f"n_loop={len(graph.loop)} n_webs={len(graph.webs)} n_offset_regions={len(regions)}")
    assert all(isinstance(r, LaminateOffsetRegion) for r in regions)
    for run in graph.loop:
        print(
            f"  loop name={run.name!r} kind={run.kind!r} "
            f"requires_eccentricity={run.requires_eccentricity} "
            f"eccentricity_m={run.eccentricity_m:.4f} T_mesh={run.thickness:.4f}"
        )
    for region in regions:
        print(
            f"  OFFSET {region.name} eccentricity_m={region.eccentricity_m:.4f} "
            f"layers={region.layers}"
        )
    strips = shell_strips_from_graph(cut)
    print(f"n_strips={len(strips)}")
    assert all(isinstance(s, ShellMidlineStrip) for s in strips)
    assert not any(s.kind == "cap" and not s.requires_eccentricity for s in strips)
    for strip in strips:
        print(
            f"  strip kind={strip.kind!r} label={strip.label!r} "
            f"n={strip.midline_twisted_pitch_axis_blade_frame.shape[0]} "
            f"T={strip.thickness_m:.4f} ecc={strip.eccentricity_m:.4f}"
        )


if __name__ == "__main__":
    main()
