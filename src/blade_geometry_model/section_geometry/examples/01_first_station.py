"""First station: NACA 2412 + layout ``2D-F``, all optional gates.

Run from the repository root::

    python src/blade_geometry_model/section_geometry/examples/01_first_station.py
    python .../01_first_station.py --plot
"""

from __future__ import annotations

import argparse
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
    plot_section,
    resolve_system_type,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Show the assembled section with plot_section",
    )
    args = parser.parse_args()

    layout = resolve_system_type("2D-F")
    airfoil = AirfoilSDF.from_naca("2412", chord=1.0)

    result = SectionGeometryRunner(
        layout=layout,
        airfoil=airfoil,
        skin_thickness=0.004,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
        twist_rad=0.0,
        outputs=("properties", "midlines", "cell_wall_loops"),
    ).execute().get_results()

    print(f"layout={layout.key}")
    print("section labels:", list(result.section.labels))
    if result.properties_map is not None:
        print("properties_map keys:", sorted(result.properties_map))
    for strip in result.midline_strips or ():
        n = strip.midline_twisted_pitch_axis_blade_frame.shape[0]
        print(f"  strip kind={strip.kind!r} label={strip.label!r} n={n}")
    for loop in result.cell_wall_loops or ():
        print(f"  cell_wall_loop id={loop.id} area_m2={loop.area_m2:.6e}")

    if args.plot:
        plot_section(result.section)


if __name__ == "__main__":
    main()
