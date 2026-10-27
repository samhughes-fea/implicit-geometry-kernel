"""Print a small ``properties_map`` table (area and centroid per label).

Run from the repository root::

    python src/blade_geometry_model/section_geometry/examples/02_read_properties.py
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


def main() -> None:
    result = (
        SectionGeometryRunner(
            layout=resolve_system_type("2D-F"),
            airfoil=AirfoilSDF.from_naca("2412", chord=1.0),
            skin_thickness=0.004,
            spar_cap_thickness=0.012,
            shear_web_thickness=0.005,
            outputs=("properties",),
        )
        .execute()
        .get_results()
    )

    props = result.properties_map
    if props is None:
        raise RuntimeError("properties_map is None; request outputs=('properties',)")

    print(f"{'label':<32} {'area_m2':>12} {'cx_m':>10} {'cy_m':>10}")
    print("-" * 66)
    for label in sorted(props):
        row = props[label]
        print(
            f"{label:<32} {row['area']:12.6e} {row['cx']:10.4f} {row['cy']:10.4f}"
        )


if __name__ == "__main__":
    main()
