"""Print the shell, figure, quadtree, and station/strip splits.

A parametric 10 m blade. NREL ``h = 1.5`` is a second block and is
skipped, with a one-line note, when the reference-blades catalogue is absent.
The figure is drawn and not saved.

    python examples/blade/16_bottleneck_splits/timings.py
"""

from __future__ import annotations

import importlib.util
import os
import runpy
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.mesh import shell_mesh  # noqa: E402
from implicit_geometry_kernel.profiling import collecting  # noqa: E402

HERE = Path(__file__).resolve().parent
REF = HERE.parents[0] / "08_reference_blades" / "reference_blades.py"
FIGURE = HERE.parents[0] / "15_shell_mesh" / "shell_mesh.py"


def parametric_blade() -> BladeGeometry:
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3, twist=8.0),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )


def _stage_seconds(reports, name: str) -> float:
    return float(sum(item.seconds for item in reports if item.stage == name))


def _shell(geo: BladeGeometry, h: float | None):
    with collecting() as handler:
        t0 = time.perf_counter()
        mesh = shell_mesh(geo) if h is None else shell_mesh(geo, h=h)
        wall = time.perf_counter() - t0
    reports = handler.take()
    return mesh, wall, _stage_seconds(reports, "shell_stations"), _stage_seconds(reports, "shell_strips")


def _figure_module():
    import matplotlib

    matplotlib.use("Agg")
    spec = importlib.util.spec_from_file_location("shell_mesh_example", FIGURE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {FIGURE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _figure(module, mesh, stem: str, title: str, text: str, window, *, save: bool, **kwargs) -> float:
    if not save:
        def _draw_only(fig, _name, folder=None, publish=False):
            fig.canvas.draw()
            module.plt.close(fig)

        module._save = _draw_only
    t0 = time.perf_counter()
    module.figure(mesh, stem, title, text, window, **kwargs)
    return time.perf_counter() - t0


def _print_shell(label: str, wall: float, figure_s: float, stations: float, strips: float) -> None:
    prefix = f"{label} " if label else ""
    print(f"{prefix}shell_mesh {wall:.3f}s  figure {figure_s:.3f}s")
    print(f"{prefix}shell_stations {stations:.3f}s  shell_strips {strips:.3f}s")


def _quadtree() -> None:
    os.environ["SECTION_GEOMETRY_PROFILE"] = "1"
    from implicit_geometry_kernel.catalogue import resolve_system_type
    from implicit_geometry_kernel.profiling import reset_stats
    from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF
    from implicit_geometry_kernel.section.assembly import build_section_view
    from implicit_geometry_kernel.section.properties.backend import quadtree_spec_for_airfoil
    from implicit_geometry_kernel.section.properties.compute import compute_section_properties

    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN-LC-0R"), twist_angle_rad=0.0)
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=7, max_cells=150000)
    reset_stats()
    with collecting() as handler:
        _props, meta = compute_section_properties(
            section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=airfoil,
        )
    phi = _stage_seconds(handler.take(), "quadtree_phi")
    profiling = meta.get("profiling") or {}
    print(f"stage_timings quadtree_phi={phi:.3f}s")
    print(f"sdf_polygon_count {int(profiling.get('sdf_polygon_count', 0))}")


def main(_argv: list[str] | None = None) -> int:
    module = _figure_module()
    geo = parametric_blade()
    mesh, wall, stations, strips = _shell(geo, None)
    rib = geo.transverse[1].z
    figure_s = _figure(
        module, mesh, "parametric", "2D-F-LC-4R, tip h/c = 5%",
        module.report("2D-F-LC-4R", mesh), (rib - 0.6, rib + 0.6),
        save=False,
    )
    _print_shell("", wall, figure_s, stations, strips)

    ns = runpy.run_path(str(REF))
    root = ns["catalogue_root"]()
    if root is None:
        print("No reference-blades catalogue; skip NREL h = 1.5.")
    else:
        _name, nrel, _meta = ns["load_windio"](ns["catalogue_yaml"](root, "nrel5mw/nrel5mw.yaml"))
        nrel_mesh, wall, stations, strips = _shell(nrel, 1.5)
        length = float(nrel.stations[-1].z)
        figure_s = _figure(
            module, nrel_mesh, "nrel5mw", "NREL 5 MW, h = 1.5 m",
            module.report("NREL 5 MW", nrel_mesh), (length - 4.0, length),
            save=False, azim=-40,
        )
        _print_shell("nrel", wall, figure_s, stations, strips)

    _quadtree()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
