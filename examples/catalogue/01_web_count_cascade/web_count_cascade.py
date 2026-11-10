"""Walk the web-count cascade and assemble each station recipe.

    python examples/catalogue/01_web_count_cascade/web_count_cascade.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.plotting import draw_cut, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.catalogue import (  # noqa: E402
    effective_layout_for_web_count,
    layout_without_spar_caps,
    resolve_system_type,
    step_down_cascade,
)

DISPLAY = 4
SKIN, CAP, WEB = 0.008, 0.03, 0.015


def main() -> None:
    parent = resolve_system_type("4D-F-LC-0R")
    ladder = step_down_cascade(parent)
    print(f"{parent.key} cascade: {ladder}")
    capless = layout_without_spar_caps(parent)
    print(f"caps run out → {capless.key} (still {capless.n_webs} webs)")
    one = effective_layout_for_web_count(parent, 1)
    print(f"one surviving web → {one.key} (D has no 1D; floor is {one.key})")

    airfoil = Airfoil.naca4("2418")
    stations = [Station(0.0, 2.0, airfoil, pitch_axis=0.3), Station(1.0, 2.0, airfoil, pitch_axis=0.3)]
    fig, axes = plt.subplots(1, len(ladder), figsize=(2.4 * len(ladder), 2.8), sharey=True)
    for ax, key in zip(axes, ladder):
        real = BladeGeometry(
            stations, key, SKIN, spar_cap_thickness=CAP, shear_web_thickness=WEB, pivot_fraction=0.3,
        ).section_at(0.5)
        cut = BladeGeometry(
            stations, key, DISPLAY * SKIN, spar_cap_thickness=DISPLAY * CAP,
            shear_web_thickness=DISPLAY * WEB, pivot_fraction=0.3,
        ).section_at(0.5)
        n_web = len(real.walls_of("web"))
        draw_cut(ax, cut)
        ax.set_title(f"{key}\n{n_web} {'web' if n_web == 1 else 'webs'}", fontsize=8)
        ax.set_xlabel("x / m")
        print(f"  assembled {key}: layout={real.layout} webs={n_web}")
    axes[0].set_ylabel("y / m")
    fig.suptitle("4D-F cascade: station recipes, parent Z/R unchanged", fontsize=10)
    note(fig, f"Display thickness {DISPLAY}×; printed web counts use the real gauges.")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, __file__, "web_count_cascade")


if __name__ == "__main__":
    main()
