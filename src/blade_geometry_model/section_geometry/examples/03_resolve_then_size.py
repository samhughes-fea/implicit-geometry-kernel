"""Resolve layout once at min skin, then rebuild with two thickness trials.

Run from the repository root::

    python src/blade_geometry_model/section_geometry/examples/03_resolve_then_size.py
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
    layout = resolve_system_type("2D-F")
    airfoil = AirfoilSDF.from_naca("2412", chord=1.0)
    t_skin_min = 0.002

    frozen = (
        SectionGeometryRunner(
            layout=layout,
            airfoil=airfoil,
            skin_thickness=t_skin_min,
            spar_cap_thickness=0.0,
            shear_web_thickness=0.0,
        )
        .resolve()
        .get_layout_resolution()
    )

    print(f"nominal={frozen.nominal_layout.key} effective={frozen.layout.key}")
    print(f"survivors_m={frozen.web_x_positions_m}")
    print(f"dropped={[(d.x_m, d.reason) for d in frozen.dropped]}")
    print(f"stepped_down={frozen.stepped_down}")

    trials = (
        (0.004, 0.012, 0.005),
        (0.006, 0.020, 0.008),
    )
    for skin_thickness, spar_cap_thickness, shear_web_thickness in trials:
        result = (
            SectionGeometryRunner(
                layout=frozen.layout,
                airfoil=airfoil,
                skin_thickness=skin_thickness,
                spar_cap_thickness=spar_cap_thickness,
                shear_web_thickness=shear_web_thickness,
                resolution=frozen,
                outputs=("properties",),
            )
            .execute()
            .get_results()
        )
        props = result.properties_map or {}
        total_area = sum(float(row["area"]) for row in props.values())
        print(
            f"trial skin={skin_thickness} spar={spar_cap_thickness} "
            f"web={shear_web_thickness} -> n_labels={len(props)} "
            f"sum_area={total_area:.6e}"
        )


if __name__ == "__main__":
    main()
