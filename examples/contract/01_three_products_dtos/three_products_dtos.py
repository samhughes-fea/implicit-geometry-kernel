"""Print ShellMidlineStrip / CellWallLoop fields and catch infeasibility.

    python examples/contract/01_three_products_dtos/three_products_dtos.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from matplotlib.lines import Line2D  # noqa: E402

from _shared.plotting import WALL_COLOURS, legend_outside, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.catalogue import effective_layout_for_web_count, resolve_system_type  # noqa: E402
from implicit_geometry_kernel.contract import (  # noqa: E402
    CellWallLoop,
    InfeasibleSectionGeometry,
    LaminateOffsetRegion,
    ShellMidlineStrip,
)
from implicit_geometry_kernel.section import (  # noqa: E402
    SectionGeometry,
    SectionOptions,
    SectionSpec,
)


def main() -> None:
    spec = SectionSpec.from_kwargs(
        resolve_system_type("2D-F-LC-0R"),
        Airfoil.naca4("2412").to_sdf(chord=1.0),
        skin_thickness=0.004,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
    )
    result = SectionGeometry(
        spec,
        SectionOptions(outputs=("properties", "midlines", "cell_wall_loops")),
    ).execute().get_results()
    strips = result.midline_strips or ()
    loops = result.cell_wall_loops or ()
    print(f"properties keys: {sorted(result.properties_map or {})}")
    print(f"n_strips={len(strips)} n_loops={len(loops)}")
    assert all(isinstance(s, ShellMidlineStrip) for s in strips)
    assert all(isinstance(c, CellWallLoop) for c in loops)
    for s in strips:
        print(f"  strip kind={s.kind!r} label={s.label!r} T={s.thickness_m:.4f} n={s.midline_twisted_pitch_axis_blade_frame.shape[0]}")
    for c in loops:
        print(f"  loop id={c.id} area_m2={c.area_m2:.6e} n={len(c.polyline_twisted_pitch_axis_blade_frame)}")
    cut = BladeGeometry(
        [Station(0.0, 1.0, Airfoil.naca4("2412"), pitch_axis=0.3),
         Station(1.0, 1.0, Airfoil.naca4("2412"), pitch_axis=0.3)],
        "2D-F-LC-0R", 0.004, spar_cap_thickness=0.012, shear_web_thickness=0.005, pivot_fraction=0.3,
    ).section_at(0.5)
    regions = cut.shell_graph().offset_regions()
    print(f"n_offset_regions={len(regions)}")
    assert all(isinstance(r, LaminateOffsetRegion) for r in regions)
    layout = resolve_system_type("2D-F-LC-0R")
    try:
        effective_layout_for_web_count(layout, 0)
    except InfeasibleSectionGeometry as exc:
        print(f"InfeasibleSectionGeometry: {exc}")

    ref = loops[0].polyline_twisted_pitch_axis_blade_frame if loops else strips[0].midline_twisted_pitch_axis_blade_frame
    swap = float(ref[:, 1].max() - ref[:, 1].min()) > float(ref[:, 0].max() - ref[:, 0].min())

    def _plot_xy(pts):
        a = pts[:, :2]
        return (a[:, 1], a[:, 0]) if swap else (a[:, 0], a[:, 1])

    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    for s in strips:
        xs, ys = _plot_xy(s.midline_twisted_pitch_axis_blade_frame)
        kind = s.kind if s.kind in WALL_COLOURS else "skin"
        ax.plot(xs, ys, color=WALL_COLOURS[kind], lw=1.4)
    for c in loops:
        xs, ys = _plot_xy(c.polyline_twisted_pitch_axis_blade_frame)
        ax.plot(xs, ys, lw=0.8, ls=":", color="#888")
    ax.set_aspect("equal")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title("contract midlines and cell loops")
    legend_outside(ax, handles=[
        Line2D([0], [0], color=WALL_COLOURS["skin"], label="skin"),
        Line2D([0], [0], color=WALL_COLOURS["web"], label="web"),
        Line2D([0], [0], color="#888", ls=":", label="cell"),
    ])
    note(fig, "SectionGeometry DTOs in the twisted pitch-axis frame: skin/web strips and one grey cell-loop set. Caps are offset regions, not midline strips.")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "three_products_dtos")


if __name__ == "__main__":
    main()
