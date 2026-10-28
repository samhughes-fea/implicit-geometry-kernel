"""
visualisation/demos/api_quickstart_demo.py
==========================================
The whole public interface, run, and what it hands back.

Every other demo reaches past the front door to show a stage of the pipeline.
This one does the opposite and shows only the front door: the twelve lines from
``section_geometry``'s README quickstart, executed verbatim, beside the section
they produce and an inventory of the result object.

The point is the ratio. Nine implicit regions, a junction-prepared centreline
topology and per-region mass properties come out of one constructor call and one
``.execute()``; there is no mesh to prepare, no boolean to repair, and no CAD
session anywhere in the loop.

The code below is deliberately *not* factored into helpers. It is meant to be
read straight off a slide and typed, so it stays flat, and it is executed here so
the published figure cannot drift away from the published code.

Writes under ``visualisation/output/api_quickstart/``:

* ``api_quickstart_section.png`` - the section, beside the result inventory.
* ``api_quickstart.json`` - the inventory, for the slide caption to quote.

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/api_quickstart_demo.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    add_output_and_presentation_args,
    apply_presentation_style,
    save_figure,
)
from blade_geometry_model.section_geometry.visualisation.palette import (
    INK,
    MUTED,
    NAVY,
    SURFACE,
    TEAL,
    strip_axes_chrome,
)


def run_quickstart():
    """The README quickstart, verbatim. Do not refactor — this is the exhibit."""
    from blade_geometry_model.section_geometry import (
        AirfoilSDF,
        SectionGeometryRunner,
        resolve_system_type,
    )

    layout = resolve_system_type("2D-F")
    airfoil = AirfoilSDF.from_naca("2412", chord=1.0)

    result = (
        SectionGeometryRunner(
            layout=layout,
            airfoil=airfoil,
            skin_thickness=0.004,
            spar_cap_thickness=0.012,
            shear_web_thickness=0.005,
            twist_rad=0.0,
            outputs=("properties", "midlines", "cell_wall_loops"),
        )
        .execute()
        .get_results()
    )
    return result


def inventory(result) -> dict:
    """What came back, in the terms the slide talks about."""
    properties = result.properties_map or {}
    strips = result.midline_strips or ()
    loops = result.cell_wall_loops or ()
    total_area = sum(float(entry["area"]) for entry in properties.values())
    return {
        "regions": len(properties),
        "midline_strips": len(strips),
        "cell_wall_loops": len(loops),
        "strip_kinds": sorted({str(strip.kind) for strip in strips}),
        "total_region_area_m2": total_area,
        "chord_m": float(result.chord_m) if result.chord_m is not None else None,
        "twist_rad": float(result.twist_rad) if result.twist_rad is not None else None,
        "property_keys": sorted(next(iter(properties.values())).keys()) if properties else [],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="api_quickstart")
    args = parser.parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    result = run_quickstart()
    facts = inventory(result)

    from blade_geometry_model.section_geometry import (
        AirfoilSDF,
        SDFGrid,
        plot_section,
    )

    # Stacked, not side by side: this sits in a half-width slide panel of about
    # 1.5:1, and a wide two-column figure letterboxes down to a thin strip there.
    figure, (axis_section, axis_result) = plt.subplots(
        2, 1, figsize=(9.0, 6.4), gridspec_kw={"height_ratios": [1.0, 1.25]}
    )

    # plot_section samples the regions, so it needs a grid. The quickstart never
    # builds one — drawing is not part of the twelve lines — so it is made here.
    grid = SDFGrid.from_airfoil(
        AirfoilSDF.from_naca("2412", chord=1.0), nx=900, ny=400
    )
    plot_section(
        result.section,
        grid,
        ax=axis_section,
        show_legend=False,
        title="",
        auto_limits=True,
    )
    axis_section.set_title(
        "…and this is what comes back", fontsize=13, color=NAVY, pad=10
    )
    strip_axes_chrome(axis_section)

    # The inventory, as a plain list. A matplotlib table would be heavier and
    # less legible at slide size than left-aligned rows.
    rows = [
        ("regions (named implicit solids)", str(facts["regions"])),
        ("midline strips (junction-prepared)", str(facts["midline_strips"])),
        ("strip kinds", ", ".join(facts["strip_kinds"])),
        ("cell-wall loops (closed, with areas)", str(facts["cell_wall_loops"])),
        ("per-region properties", ", ".join(facts["property_keys"][:6])),
        ("chord", f"{facts['chord_m']:.3f} m"),
        ("lines of code above", "12"),
    ]

    axis_result.axis("off")
    axis_result.set_title("result: SectionGeometryResult", fontsize=13, color=NAVY, pad=10)
    top = 0.93
    step = 0.135
    for index, (label, value) in enumerate(rows):
        y = top - index * step
        axis_result.text(0.0, y, label, fontsize=12, color=MUTED, va="center")
        axis_result.text(
            1.0,
            y,
            value,
            fontsize=13.5,
            color=TEAL if index < len(rows) - 1 else INK,
            fontweight="bold",
            ha="right",
            va="center",
        )
        axis_result.plot(
            [0.0, 1.0],
            [y - 0.5 * step, y - 0.5 * step],
            color="#e3e7ec",
            linewidth=0.8,
            clip_on=False,
        )
    axis_result.set_xlim(0.0, 1.0)
    axis_result.set_ylim(0.0, 1.0)
    axis_result.set_facecolor(SURFACE)

    figure.suptitle(
        "No mesh to prepare, no boolean to repair, no CAD session in the loop",
        fontsize=13.5,
    )
    figure.tight_layout()
    out_png = out_dir / "api_quickstart_section.png"
    save_figure(figure, out_png, presentation=presentation)
    plt.close(figure)
    print(f"Wrote {out_png}")

    out_json = out_dir / "api_quickstart.json"
    out_json.write_text(json.dumps(facts, indent=2), encoding="utf-8")
    print(f"Wrote {out_json}")
    for label, value in rows:
        print(f"  {label:<38} {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
