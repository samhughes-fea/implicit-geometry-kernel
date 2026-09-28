"""Handoff midlines to ``ShellStripMeshingInputs`` (stop before a full mesh study).

Optionally builds a coarse strip mesh when ``--mesh`` is passed.

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/examples/05_handoff_to_mesh.py
    python .../05_handoff_to_mesh.py --mesh
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

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    SectionGeometryRunner,
    resolve_system_type,
)
from blade_structure_utilities.mesh_utilities.section_mesh import (
    shell_mesh_inputs_from_geometry_result,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mesh",
        action="store_true",
        help="Also build a coarse SectionStripMesh (n_elements_per_panel=4)",
    )
    args = parser.parse_args()

    layout_key = "2D-F"
    result = (
        SectionGeometryRunner(
            layout=resolve_system_type(layout_key),
            airfoil=AirfoilSDF.from_naca("2412", chord=1.0),
            skin_thickness=0.004,
            spar_cap_thickness=0.012,
            shear_web_thickness=0.005,
            outputs=("midlines",),
        )
        .execute()
        .get_results()
    )

    shell = shell_mesh_inputs_from_geometry_result(result, layout_key=layout_key)
    print(f"layout_key={shell.layout_key}")
    print(f"chord_m={shell.chord_m} twist_rad={shell.twist_rad}")
    print(f"n_midlines={len(shell.midlines)}")
    print(f"LE={shell.leading_edge_twisted_pitch_axis_blade_frame}")
    print(f"TE={shell.trailing_edge_twisted_pitch_axis_blade_frame}")
    for m in shell.midlines:
        n = m.midline_twisted_pitch_axis_blade_frame.shape[0]
        print(f"  {m.kind!r} {m.label!r} n={n} t={m.thickness_m}")

    if args.mesh:
        from blade_structure_utilities.mesh_utilities.section_mesh import (
            build_section_strip_mesh,
        )

        mesh = build_section_strip_mesh(
            shell, n_elements_per_panel=4, endpoint_tol=1e-6
        )
        print(
            f"strip_mesh panels={len(mesh.panel_meshes)} "
            f"n_total_elements={mesh.n_total_elements}"
        )


if __name__ == "__main__":
    main()
