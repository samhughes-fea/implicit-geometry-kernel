"""Time the shell meshers on a generic blade and write the numbers as JSON.

    python scripts/bench_shell_mesh.py                  # default sizes, prints JSON
    python scripts/bench_shell_mesh.py --h 0.2 0.1 --out bench.json

Each row: mesher (``ruled`` or ``partition``), element size ``h``, wall time to
mesh, time to write binary ``.vtu`` and legacy ``.vtk``, node and element counts
and the mesh digest. The blade is a tapered, twisted NACA 0018 with two webs and
two ribs, built from public inputs only. Numbers depend on the machine: compare
rows from one run, not across machines.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from pathlib import Path


def blade():
    from implicit_geometry_kernel import Airfoil, BladeGeometry, Station

    foil = Airfoil.naca4("0018")
    stations = [
        Station(0.0, 1.0, foil, twist=5.0, pitch_axis=0.3),
        Station(4.0, 0.7, foil, twist=0.0, pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-2R", 0.004, spar_cap_thickness=0.008,
                         shear_web_thickness=0.004, pivot_fraction=0.3)


def run(sizes, methods=("ruled", "partition")) -> list[dict]:
    from implicit_geometry_kernel.mesh import shell_mesh, write_vtk, write_vtu

    rows = []
    for h in sizes:
        for method in methods:
            geo = blade()  # fresh geometry: no cut cache carried between rows
            t0 = time.perf_counter()
            mesh = shell_mesh(geo, h=float(h), method=method)
            t_mesh = time.perf_counter() - t0
            with tempfile.TemporaryDirectory() as tmp:
                t0 = time.perf_counter()
                write_vtu(mesh, Path(tmp) / "m.vtu")
                t_vtu = time.perf_counter() - t0
                t0 = time.perf_counter()
                write_vtk(mesh, Path(tmp) / "m.vtk")
                t_vtk = time.perf_counter() - t0
            rows.append({
                "mesher": method, "h": float(h), "mesh_s": round(t_mesh, 3), "vtu_s": round(t_vtu, 4),
                "vtk_s": round(t_vtk, 4), "nodes": len(mesh.nodes), "elements": len(mesh.elements),
                "digest": mesh.digest(),
            })
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--h", type=float, nargs="+", default=[0.2, 0.1])
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)
    text = json.dumps(run(args.h), indent=2)
    if args.out is None:
        sys.stdout.write(text + "\n")
    else:
        args.out.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
