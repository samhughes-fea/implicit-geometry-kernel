"""
Registry sweep: build section geometry for every ``SYSTEM_TYPE_KEYS`` layout and
save SDF subcomponent figures + taxonomy coverage JSON.

Owns the **geometry** half of the former FE ``repro_section_centerline_mesh_v2``
script (40-key airfoil/layout/section generation). Mesh strip figures live under
``section_mesh/visualisation/demos/``.

Outputs (under ``visualisation/output/section_centerline_geometry_taxonomy/``):

* ``section_geom_<layout>_i000.png`` — SDF subcomponents (skin/web/cap)
* ``taxonomy_coverage.json`` — rendered vs full registry

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/section_centerline_geometry_taxonomy.py
    python .../section_centerline_geometry_taxonomy.py 0A 0B 2D-F
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    demo_output_dir,
    save_figure,
)


def _all_layout_keys() -> tuple[str, ...]:
    from blade_geometry_model.section_geometry.layout.system_layout import (
        SYSTEM_TYPE_KEYS,
    )

    return SYSTEM_TYPE_KEYS


def _expected_registry_keys() -> set[str]:
    return {"0A", "0B"} | {
        f"{x}{y}-{z}"
        for x in range(1, 6)
        for y in ("A", "B", "C", "D")
        for z in ("CN", "F")
        if not (x == 1 and y == "D")
    }


def _save_subcomponent_figure(
    *,
    outfile: Path,
    section: Any,
    airfoil: Any,
    layout: Any,
    layout_key: str,
    twist_rad: float,
    twist_deg: float,
    strips_geo: tuple[Any, ...],
    station_tag: str = "NACA63-415",
    dpi: int = 160,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import numpy as np

    from blade_geometry_model.coordinate_frames import (
        naca_to_twisted_pitch_axis_blade_frame,
        section_plot_axes,
    )
    from blade_geometry_model.section_geometry.visualisation.plot_section import (
        apply_structural_section_axes,
        fit_square_section_axes,
    )
    from blade_geometry_model.section_geometry.visualisation.plot_subcomponents import (
        draw_section_sdf_panel,
        taxonomy_plot_caption,
    )

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#0f1117")

    draw_section_sdf_panel(ax, section, airfoil, twist_rad)

    # LE/TE in canonical B (same handoff as midlines), not bare airfoil.rotate.
    chord_m = float(airfoil.chord)
    le_te_naca = np.asarray([airfoil.leading_edge, airfoil.trailing_edge], dtype=float)
    le_te_b = naca_to_twisted_pitch_axis_blade_frame(
        le_te_naca, float(twist_rad), chord_m
    )
    h, v = section_plot_axes(le_te_b)
    ax.scatter(
        h,
        v,
        color="#ffffff",
        s=60,
        marker="x",
        zorder=30,
        label="LE / TE",
    )

    apply_structural_section_axes(ax)
    fit_square_section_axes(ax)
    ax.tick_params(colors="#aaaaaa")
    ax.grid(True, color="#2a2a3a", lw=0.4, alpha=0.5)
    ax.xaxis.label.set_color("#e0e0e0")
    ax.yaxis.label.set_color("#e0e0e0")

    tax = taxonomy_plot_caption(layout_key, layout)
    n_sub = len(getattr(section, "labels", []) or [])
    n_strips = len(strips_geo)
    tag_txt = f" · {station_tag}" if station_tag else ""
    fig.suptitle(
        f"{layout_key}{tag_txt} — geometry: Subcomponents (SDF)\n"
        f"{tax}\n"
        f"twist = {twist_deg:.2f}° · n_subcomponents={n_sub} · n_midline_strips={n_strips}",
        color="#e0e0e0",
        fontsize=11,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90))
    save_figure(fig, outfile, dpi=dpi, facecolor="#0f1117")
    plt.close(fig)
    return outfile


def _run_single_layout(
    *,
    layout_key: str,
    airfoil: Any,
    twist_rad: float,
    twist_deg: float,
    out_dir: Path,
) -> None:
    from blade_geometry_model.section_geometry.api.factory import (
        build_section_view,
        resolve_system_type,
    )
    from blade_geometry_model.section_geometry.midlines.shell_midline_export import (
        build_shell_midline_strips,
    )

    layout = resolve_system_type(layout_key)
    section = build_section_view(airfoil, layout, twist_angle_rad=twist_rad)
    strips_geo = build_shell_midline_strips(
        section, twist_rad=twist_rad, n_web_samples=20, n_cap_samples=80
    )
    assert strips_geo, f"[{layout_key}] expected midline strips from geometry export"

    fig_path = out_dir / f"section_geom_{layout_key}_i000.png"
    _save_subcomponent_figure(
        outfile=fig_path,
        section=section,
        airfoil=airfoil,
        layout=layout,
        layout_key=layout_key,
        twist_rad=twist_rad,
        twist_deg=twist_deg,
        strips_geo=strips_geo,
    )
    print(f"[section_centerline_geometry] wrote {fig_path} ({len(strips_geo)} strips)")


def main(argv: list[str] | None = None) -> None:
    import numpy as np

    from blade_geometry_model.section_geometry.api.factory import (
        resolve_system_type,
    )
    from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF

    argv = list(sys.argv[1:] if argv is None else argv)

    # Same station fixture as the former FE centerline repro (run007 / i000).
    chord_m = 1.655959896
    twist_deg = 34.110405083
    airfoil = AirfoilSDF.from_naca_series(
        6, 63.0, 4.0, 15.0, n_points=200, chord=chord_m, closed_te=True
    )
    twist_rad = float(np.deg2rad(twist_deg))

    out_dir = (
        Path(__file__).resolve().parents[1] / "output" / "section_centerline_geometry_taxonomy"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    if argv:
        layout_keys = [resolve_system_type(k.strip()).key for k in argv if k.strip()]
    else:
        layout_keys = list(_all_layout_keys())

    if not layout_keys:
        raise SystemExit("No layout keys to render.")

    rendered_png_names = {f"section_geom_{k}_i000.png" for k in layout_keys}
    stale_pngs = sorted(
        p.name for p in out_dir.glob("section_geom_*_i000.png") if p.name not in rendered_png_names
    )
    for name in stale_pngs:
        try:
            (out_dir / name).unlink()
        except OSError:
            pass

    expected_keys = _expected_registry_keys()
    rendered_set = set(layout_keys)
    coverage_manifest = {
        "rendered_keys": list(layout_keys),
        "missing_vs_full_registry": sorted(expected_keys - rendered_set),
        "extra_vs_full_registry": sorted(rendered_set - expected_keys),
        "excluded_invalid_keys": ["1D-CN", "1D-F"],
        "stale_pngs_removed": stale_pngs,
    }
    with open(out_dir / "taxonomy_coverage.json", "w", encoding="utf-8") as f:
        json.dump(coverage_manifest, f, indent=2)

    print(
        f"[section_centerline_geometry] rendering {len(layout_keys)} layout(s): "
        f"{', '.join(layout_keys)}"
    )
    for layout_key in layout_keys:
        _run_single_layout(
            layout_key=layout_key,
            airfoil=airfoil,
            twist_rad=twist_rad,
            twist_deg=twist_deg,
            out_dir=out_dir,
        )


if __name__ == "__main__":
    main()
