"""Assemble a thin-walled box from structure contributions and span laws.

    python examples/structure/03_assemble_box_beam/assemble_box_beam.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from _shared.plotting import WALL_COLOURS, legend_outside, meta_cmap, note, plt, save  # noqa: E402
from implicit_geometry_kernel.sdf.nodes.tagged import PolygonSDF  # noqa: E402
from implicit_geometry_kernel.structure import (  # noqa: E402
    BranchContribution,
    BuildContext,
    ConstantLaw,
    TableLaw,
    WallContribution,
    assemble,
)

W, H, L = 1.0, 0.2, 6.0
SKIN, CAP, WEB = 0.01, 0.03, 0.008


def rect(x0, x1, y0, y1) -> PolygonSDF:
    return PolygonSDF(np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]]))


def main() -> None:
    ctx = BuildContext(oml=rect(-W / 2, W / 2, -H / 2, H / 2), span=(0.0, L))
    skin = WallContribution(ctx.table.add("skin", "skin"), TableLaw([0.0, L], [SKIN, 0.6 * SKIN]))
    cap = WallContribution(ctx.table.add("cap_upper", "cap"), ConstantLaw(CAP),
                           footprint=rect(-0.2, 0.2, 0.0, 1.0))
    web = BranchContribution(
        ctx.table.add("web", "web"),
        rect(-1e-9, 1e-9, -1.0, 1.0),
        rect(-WEB / 2, WEB / 2, -1.0, 1.0),
        ConstantLaw(WEB),
    )
    structure = assemble(ctx, [skin, cap, web])
    z = 2.0
    pts = np.array([[0.4, H / 2, z], [0.0, H / 2, z], [0.0, 0.0, z]])
    t = structure.thickness(pts)
    print(f"T(skin, cap, web-plane) = {t}")
    for x, y, label in ((0.4, H / 2 - 0.5 * SKIN, "skin"), (0.0, H / 2 - 0.5 * CAP, "cap_upper"), (0.0, 0.0, "web")):
        s = structure.material.evaluate(np.array([[x, y, z]]))
        name = structure.table[int(s.meta[0])].name if s.phi[0] < 0 else None
        print(f"  ({x:.2f}, {y:.3f}): {name} (want {label})")

    names = {spec.id: spec.name for spec in structure.table}
    colours = [WALL_COLOURS.get(structure.table[i].kind, "#888") for i in range(len(structure.table))]
    cmap, norm = meta_cmap(len(structure.table), colours)

    xs = np.linspace(-W / 2 - 0.05, W / 2 + 0.05, 160)
    ys = np.linspace(-H / 2 - 0.05, H / 2 + 0.05, 80)
    X, Y = np.meshgrid(xs, ys)
    q = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, z)])
    sample = structure.material.evaluate(q)
    owner = np.where(sample.phi < 0.0, sample.meta, np.nan).reshape(X.shape)

    fig, ax = plt.subplots(figsize=(7.4, 2.9))
    ax.pcolormesh(X, Y, owner, cmap=cmap, norm=norm, shading="nearest")
    ax.plot([-W / 2, W / 2, W / 2, -W / 2, -W / 2],
            [-H / 2, -H / 2, H / 2, H / 2, -H / 2], color="k", lw=0.6, ls="--")
    ax.set_aspect("equal")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title("assembled material at z = 2 m")
    legend_outside(ax, handles=[Patch(facecolor=colours[i], label=names[i]) for i in sorted(names)])
    note(fig, "skin navy, cap_upper red, web teal. Outside the walls is empty (white).")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "assemble_box_beam")


if __name__ == "__main__":
    main()
