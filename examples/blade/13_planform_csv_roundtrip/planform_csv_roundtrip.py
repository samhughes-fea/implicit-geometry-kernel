"""read_planform_csv → BladeGeometry → write_windio.

    python examples/blade/13_planform_csv_roundtrip/planform_csv_roundtrip.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.plotting import draw_shell, note, plt, save  # noqa: E402
from implicit_geometry_kernel import read_planform_csv, write_windio  # noqa: E402

HERE = Path(__file__).resolve().parent


def main() -> None:
    geo = read_planform_csv(
        HERE / "planform.csv",
        layout="2D-F-LC-0R",
        skin_thickness=0.004,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.004,
    )
    doc = write_windio(geo)
    bem = doc["outer_shape_bem"]
    print(f"stations={len(geo.stations)} span={geo.stations[-1].z - geo.stations[0].z:.1f} m")
    print("windIO chords:", bem["chord"]["values"])
    print("windIO labels:", bem["airfoil_position"]["labels"])
    zs = (0.0, 5.0, 10.0)
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.7), sharey=True)
    for ax, z in zip(axes, zs):
        draw_shell(ax, geo.section_at(z))
        ax.set_title(f"z = {z:g} m", fontsize=8)
        ax.set_xlabel("x / m")
    axes[0].set_ylabel("y / m")
    fig.suptitle("planform CSV loft", fontsize=10)
    note(fig, "Shared y-limits: the tip is small because the CSV chord tapers.")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, __file__, "planform_csv_roundtrip")


if __name__ == "__main__":
    main()
