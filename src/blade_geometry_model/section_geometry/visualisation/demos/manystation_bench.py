"""
visualisation/demos/manystation_bench.py
=============================
Many-station throughput benchmark for the section_geometry property backend.

Times the per-station hot path (section build + property integration,
and optionally plot + JSON export) across ``--stations`` stations.

Run from the repository root:
    python src/blade_geometry_model/section_geometry/visualisation/demos/manystation_bench.py --stations 25
    python src/blade_geometry_model/section_geometry/visualisation/demos/manystation_bench.py --stations 10 --backend hybrid
    python src/blade_geometry_model/section_geometry/visualisation/demos/manystation_bench.py --stations 10 --compare
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import tempfile
import time
from pathlib import Path

# Put the repo root on sys.path before the first package import, as every other
# demo does. Without it this script only ran from inside ``src/``, and the
# invocation in its own docstring raised ModuleNotFoundError.
_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)

# Enable phi-eval profiling before engine imports.
os.environ.setdefault("SECTION_GEOMETRY_PROFILE", "1")

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
if _erp.exists():
    _spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
    _mod = importlib.util.module_from_spec(_spec)
    assert _spec.loader is not None
    _spec.loader.exec_module(_mod)
    _mod.ensure_repo_path()

import numpy as np

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.sdf.cache import SectionEvalCache
from blade_geometry_model.section_geometry.sdf.profiling import get_stats, reset_stats


def _station_params(n: int):
    """Spread stations over a few NACA codes and twist angles."""
    codes = ["2412", "2410", "4412", "0012"]
    for i in range(n):
        code = codes[i % len(codes)]
        twist = np.deg2rad(15.0 * (i / max(1, n - 1)))  # 0 .. 15 deg
        yield code, float(twist)


def _run_bench(args, backend: str) -> dict:
    sys_type = resolve_system_type(args.system)
    tmp = Path(tempfile.mkdtemp(prefix=f"manystation_bench_{backend}_"))
    stage_s = {"build": 0.0, "properties": 0.0, "plot": 0.0, "export": 0.0}

    reset_stats()
    t0 = time.perf_counter()
    last_meta = {}
    for i, (code, twist) in enumerate(_station_params(args.stations)):
        t = time.perf_counter()
        airfoil = AirfoilSDF.from_naca(code, n_points=args.n_points, chord=1.0)
        section = build_section_view(airfoil, sys_type, twist_angle_rad=twist)
        stage_s["build"] += time.perf_counter() - t

        spec = quadtree_spec_for_airfoil(
            airfoil, chord=1.0, max_depth=args.depth, max_cells=150000
        )
        eval_cache = SectionEvalCache()

        t = time.perf_counter()
        _props, _meta, mesh = compute_section_properties(
            section,
            backend=backend,
            quadtree_spec=spec,
            airfoil_for_spec=airfoil,
            with_mesh=True,
        )
        stage_s["properties"] += time.perf_counter() - t
        last_meta = _meta

        if args.plot and mesh is not None:
            from blade_geometry_model.section_geometry.visualisation.plot import (
                plot_section_from_mesh,
                plot_section_hybrid,
            )
            import matplotlib

            matplotlib.use("Agg", force=False)
            import matplotlib.pyplot as plt

            t = time.perf_counter()
            if backend == "hybrid":
                fig, _ = plot_section_hybrid(section, mesh)
            else:
                fig, _ = plot_section_from_mesh(mesh)
            save_figure(fig, tmp / f"station_{i:03d}.png")
            plt.close(fig)
            stage_s["plot"] += time.perf_counter() - t

        if args.export:
            from blade_geometry_model.section_geometry.props_emit.export import (
                SectionPropertiesReport,
                export_section_json,
            )
            from blade_geometry_model.section_geometry.sdf.grid import SDFGrid

            grid = SDFGrid.from_bbox(spec.x_min, spec.x_max, spec.y_min, spec.y_max, 384, 180)
            report = SectionPropertiesReport.from_precomputed(
                section, _props, _meta, mesh
            )
            t = time.perf_counter()
            export_section_json(
                section,
                grid,
                {},
                str(tmp / f"station_{i:03d}.json"),
                include_component_zero_contours=True,
                include_geometry_detail=True,
                eval_cache=eval_cache,
                properties_report=report,
            )
            stage_s["export"] += time.perf_counter() - t

    wall = time.perf_counter() - t0
    n = args.stations
    stats = get_stats()
    return {
        "backend": backend,
        "wall": wall,
        "n": n,
        "stage_s": stage_s,
        "stats": stats.to_dict() if stats is not None else {},
        "meta": last_meta,
        "tmp": tmp,
        "ms_per_station": wall / n * 1000,
        "properties_ms": stage_s["properties"] / n * 1000,
    }


def _print_result(r: dict) -> None:
    n = r["n"]
    print(f"\nMany-station bench: {n} stations, backend={r['backend']}")
    print(f"  total wall        : {r['wall']:8.3f} s   ({r['ms_per_station']:7.1f} ms/station)")
    for stage, v in r["stage_s"].items():
        if v:
            print(f"  {stage + '_s':<18}: {v:8.3f} s   ({v / n * 1000:7.1f} ms/station)")
    d = r["stats"]
    if d:
        print(f"  phi_eval_count    : {d.get('phi_eval_count', 0):>12,}")
        print(f"  sdf_polygon_count : {d.get('sdf_polygon_count', 0):>12,}")
    meta = r["meta"]
    if meta.get("constructive_s") is not None:
        print(f"  constructive_s    : {meta['constructive_s']:.4f} (last station)")
        print(f"  residual_quadtree_s: {meta.get('residual_quadtree_s', 0):.4f} (last station)")
        print(f"  residual_labels   : {meta.get('residual_labels')}")
    print(f"  artifacts in      : {r['tmp']}")


def _plot_throughput(
    results: list[dict],
    path: Path,
    *,
    stations_per_blade: int,
    designs: int,
    presentation: bool,
) -> None:
    """Per-station cost, and what it extrapolates to at blade and sweep scale.

    The right-hand panel is an extrapolation, not a measurement, and is labelled
    as one: it multiplies the measured per-station cost by a stated station
    count. That is the honest way to answer "can this run inside a co-design
    loop?", which is a question about blades and sweeps, not about one station.
    """
    import matplotlib

    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    from blade_geometry_model.section_geometry.visualisation.palette import (
        INK,
        MUTED,
        NAVY,
        TEAL,
    )

    figure, (axis_stage, axis_scale) = plt.subplots(1, 2, figsize=(14.0, 5.2))

    backends = [str(r["backend"]) for r in results]
    positions = np.arange(len(results))
    build_ms = [1000.0 * r["stage_s"]["build"] / r["n"] for r in results]
    properties_ms = [1000.0 * r["stage_s"]["properties"] / r["n"] for r in results]

    axis_stage.bar(positions, build_ms, width=0.55, color=NAVY, label="build")
    axis_stage.bar(
        positions, properties_ms, width=0.55, bottom=build_ms, color=TEAL, label="properties"
    )
    for position, result in zip(positions, results):
        axis_stage.text(
            position,
            float(result["ms_per_station"]) * 1.03,
            f"{float(result['ms_per_station']):.0f} ms",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            color=INK,
        )
    axis_stage.set_xticks(positions)
    axis_stage.set_xticklabels(backends)
    axis_stage.set_ylabel("wall time per station [ms]")
    axis_stage.set_title("Where a station's time goes")
    axis_stage.grid(True, axis="y", alpha=0.25)
    axis_stage.legend(fontsize=9)

    fastest = min(results, key=lambda r: float(r["ms_per_station"]))
    seconds_per_station = float(fastest["ms_per_station"]) / 1000.0
    tiers = [
        ("one station", 1),
        (f"one blade\n({stations_per_blade} stations)", stations_per_blade),
        (
            f"a {designs}-design sweep\n({designs * stations_per_blade:,} sections)",
            designs * stations_per_blade,
        ),
    ]
    tier_positions = np.arange(len(tiers))
    tier_seconds = [seconds_per_station * count for _, count in tiers]
    axis_scale.bar(
        tier_positions, tier_seconds, width=0.55, color=NAVY
    )
    for position, seconds in zip(tier_positions, tier_seconds):
        text = f"{seconds * 1000:.0f} ms" if seconds < 1.0 else f"{seconds:.1f} s"
        axis_scale.text(
            position,
            seconds * 1.12,
            text,
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            color=INK,
        )
    axis_scale.set_xticks(tier_positions)
    axis_scale.set_xticklabels([label for label, _ in tiers], fontsize=10)
    axis_scale.set_yscale("log")
    axis_scale.set_ylabel("wall time [s]")
    axis_scale.set_title(
        f"Extrapolated from {fastest['backend']} at "
        f"{float(fastest['ms_per_station']):.0f} ms/station, single core"
    )
    axis_scale.grid(True, axis="y", alpha=0.25)
    # Under the axis rather than inside it: floated in the panel this caption sat
    # across the "one station" bar.
    axis_scale.set_xlabel(
        "extrapolation: measured per-station cost × station count",
        fontsize=8.5,
        style="italic",
        color=MUTED,
    )

    figure.suptitle(
        "Geometry at co-design speed: a whole blade of sections in "
        f"{seconds_per_station * stations_per_blade:.1f} s",
        fontsize=13,
    )
    figure.tight_layout()
    save_figure(figure, path, presentation=presentation)
    plt.close(figure)
    print(f"Wrote {path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stations", type=int, default=20)
    ap.add_argument("--depth", type=int, default=7)
    ap.add_argument("--system", default="2D-CN")
    ap.add_argument("--n-points", type=int, default=120)
    ap.add_argument("--backend", default="hybrid", choices=("hybrid", "quadtree", "grid"))
    ap.add_argument(
        "--compare",
        action="store_true",
        help="run both hybrid and quadtree and print speedup",
    )
    ap.add_argument("--plot", action="store_true", help="also render a PNG per station")
    ap.add_argument("--export", action="store_true", help="also write a section JSON per station")
    ap.add_argument(
        "--figure",
        action="store_true",
        help="write the throughput figure (implies --compare).",
    )
    ap.add_argument(
        "--stations-per-blade",
        type=int,
        default=40,
        help="Station count used for the extrapolated blade tier in --figure.",
    )
    ap.add_argument(
        "--designs",
        type=int,
        default=100,
        help="Design count used for the extrapolated sweep tier in --figure.",
    )
    ap.add_argument(
        "--presentation",
        action="store_true",
        help="Larger fonts and higher DPI for slides/decks.",
    )
    args = ap.parse_args()
    apply_presentation_style(bool(args.presentation))

    if args.compare or args.figure:
        r_qt = _run_bench(args, "quadtree")
        _print_result(r_qt)
        r_hy = _run_bench(args, "hybrid")
        _print_result(r_hy)
        speedup = r_qt["properties_ms"] / max(r_hy["properties_ms"], 1e-9)
        print(f"\n  properties speedup (quadtree/hybrid): {speedup:.2f}x")
        print(f"  hybrid properties_ms/station: {r_hy['properties_ms']:.1f}")
        if args.figure:
            _plot_throughput(
                [r_hy, r_qt],
                demo_output_dir("manystation_bench") / "manystation_throughput.png",
                stations_per_blade=int(args.stations_per_blade),
                designs=int(args.designs),
                presentation=bool(args.presentation),
            )
        return

    r = _run_bench(args, args.backend)
    _print_result(r)


if __name__ == "__main__":
    main()
