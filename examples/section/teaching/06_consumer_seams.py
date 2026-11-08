"""Consume Step 6 seams the way a laminate or section-model repo should.

Gmsh jobs, ply ABD/PCOMP and solvers stay downstream. This kernel only hands
over metres: ``ThicknessLaw``, ``named_edges``, ``size_fields``, and
``merge_cap_and_skin``. Clip with ``apply_clip``. Pin
``implicit-geometry-kernel>=0.8.2`` if you read indexed C-family pads or mapped
loft faces.

Run from the repository root::

    python examples/section/teaching/06_consumer_seams.py
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

import numpy as np

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.stations import PlyDrop
from implicit_geometry_kernel.structure.laws import ThicknessLaw


def main() -> None:
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    drop = PlyDrop("d1", 0.2, 0.4, 1.0, 5.0, 0.002, taper=0.0)
    stacked = BladeGeometry(
        stations,
        "2D-F-LC-0R",
        0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48),
        ply_drops=(drop,),
        pivot_fraction=0.3,
    )
    merged = BladeGeometry(
        stations,
        "2D-F-LC-0R",
        0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48),
        merge_cap_and_skin=True,
        pivot_fraction=0.3,
    )
    law: ThicknessLaw = stacked.thickness_law()
    z, s = 3.0, 0.3
    print(f"thickness_law(z={z}, s={s}) = {float(law(z, s)):.6f} m")
    print(f"forced_stations = {stacked.forced_stations()}")
    edges = stacked.named_edges(4.0)
    print(f"n_named_edges @ z=4 = {len(edges)}")
    for e in edges:
        if "d1:" in f"{e.name_a}{e.name_b}":
            print(f"  ply-drop edge {e.name_a!r} {e.name_b!r}")
    cut = stacked.section_at(4.0)
    mid = cut.shell_graph().loop[0].points[len(cut.shell_graph().loop[0].points) // 2]
    fields = stacked.size_fields(mid.reshape(1, 3))
    print("size_fields keys:", sorted(fields))
    print(f"stacked walls = {len(stacked.section_at(4.0).walls)}")
    print(f"merged walls  = {len(merged.section_at(4.0).walls)} (graph runs; cap thickness stacked)")
    cap = next(w for w in merged.section_at(4.0).walls if w.kind == "cap")
    print(f"merged cap thickness = {cap.thickness:.6f} m")
    np.testing.assert_allclose(law(z, s), stacked.thickness_at(z, s))


if __name__ == "__main__":
    main()
