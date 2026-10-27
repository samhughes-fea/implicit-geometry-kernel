"""
visualisation/demos/constructive_vs_quadtree_spike.py
=========================================
WP1 spike: constructive wall-solid polygons + shoelace vs full-label quadtree.

Reports emit/shoelace timings and per-label area errors for 2D-CN NACA 2412.
Naive buffer is not used — emitters are recipe-faithful (WP2).

Run from repository root:
    python src/blade_geometry_model/section_geometry/visualisation/demos/constructive_vs_quadtree_spike.py
"""

from __future__ import annotations

import importlib.util
import time
from pathlib import Path

import numpy as np

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.constructive import (
    emit_constructive_polygons,
    is_constructive_label,
    section_properties_from_polygons,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF


def main() -> None:
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=7, max_cells=150000)

    t0 = time.perf_counter()
    polys = emit_constructive_polygons(section)
    emit_s = time.perf_counter() - t0
    t1 = time.perf_counter()
    _props = section_properties_from_polygons(polys)
    shoelace_s = time.perf_counter() - t1
    del _props

    t2 = time.perf_counter()
    qt, qt_meta = compute_section_properties(
        section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    qt_s = time.perf_counter() - t2

    t3 = time.perf_counter()
    hy, hy_meta = compute_section_properties(
        section, backend="hybrid", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    hy_s = time.perf_counter() - t3

    print("=== constructive vs quadtree spike (2D-CN NACA 2412) ===")
    print(f"emit_s={emit_s:.4f}  shoelace_s={shoelace_s:.6f}  constructive_total={emit_s + shoelace_s:.4f}")
    print(f"quadtree_s={qt_s:.4f}  hybrid_s={hy_s:.4f}  speedup={qt_s / max(hy_s, 1e-9):.2f}x")
    print(f"hybrid constructive_s={hy_meta.get('constructive_s')} residual_s={hy_meta.get('residual_quadtree_s')}")
    print(f"residual_labels={hy_meta.get('residual_labels')}")
    print()
    print(f"{'label':<18} {'qt_A':>10} {'hy_A':>10} {'rel_err':>10}")
    for label in section.labels:
        key = str(getattr(label, "value", label))
        if not is_constructive_label(key):
            continue
        a_qt = float(qt[key]["area"])
        a_hy = float(hy[key]["area"])
        err = abs(a_hy - a_qt) / max(a_qt, 1e-15)
        print(f"{key:<18} {a_qt:10.5f} {a_hy:10.5f} {err:10.3e}")

    print("\nCores are constructive; LE/TE inserts remain on residual-phi when enabled.")


if __name__ == "__main__":
    main()
