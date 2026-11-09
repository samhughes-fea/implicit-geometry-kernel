"""Compare ring_sums_S, quadtree quadrature, and an exact disk.

    python examples/integration/01_quadtree_vs_shoelace/quadtree_vs_shoelace.py
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

from _shared.plotting import SERIES, note, plt, save  # noqa: E402
from implicit_geometry_kernel.integration import (  # noqa: E402
    QuadtreeSpec,
    build_shared_quadtree,
    ring_sums_S,
    section_properties_from_quadtree,
)
from implicit_geometry_kernel.sdf import sdf_circle  # noqa: E402

R = 0.30
EXACT = math.pi * R * R


def shoelace_area(n: int) -> float:
    t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    ring = np.column_stack([R * np.cos(t), R * np.sin(t)])
    cross = ring_sums_S(ring)[0]
    return 0.5 * abs(cross)


def as_sampler(fn):
    return lambda pts: fn(np.asarray(pts, float)[:, 0], np.asarray(pts, float)[:, 1])


def quadtree_area(h: float) -> float:
    spec = QuadtreeSpec(-0.5, 0.5, -0.5, 0.5, max_depth=14, min_cell_size=h, refine_band=2 * h,
                        max_cells=500_000, root_cell_size=0.1, refine_probes="full")
    disk = lambda x, y: sdf_circle(x, y, 0.0, 0.0, R)  # noqa: E731
    mesh = build_shared_quadtree(spec, {"disk": as_sampler(disk)})
    props = section_properties_from_quadtree(mesh, mesh.phi_by_label, "disk", sampler=as_sampler(disk))
    return float(props["area"])


def main() -> None:
    ns = (16, 32, 64, 128)
    hs = (0.08, 0.04, 0.02)
    shoe = [shoelace_area(n) for n in ns]
    qt = [quadtree_area(h) for h in hs]
    print(f"exact area {EXACT:.6f} m^2")
    for n, a in zip(ns, shoe):
        print(f"  shoelace n={n:3d}: {a:.6f}  rel err {abs(a - EXACT) / EXACT:.2e}")
    for h, a in zip(hs, qt):
        print(f"  quadtree h={h:.2f}: {a:.6f}  rel err {abs(a - EXACT) / EXACT:.2e}")

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4), sharey=True)
    axes[0].semilogy(ns, [abs(a - EXACT) / EXACT for a in shoe], "o-", color=SERIES[0])
    axes[0].set_xlabel("n vertices on the ring")
    axes[0].set_ylabel("relative area error")
    axes[0].set_title("shoelace (Green) on a sampled circle")
    axes[1].semilogy(np.array(hs) * 1e3, [abs(a - EXACT) / EXACT for a in qt], "s-", color=SERIES[1])
    axes[1].set_xlabel("min cell size h / mm (left = coarser)")
    axes[1].invert_xaxis()
    axes[1].set_title("quadtree quadrature of the disk SDF")
    fig.suptitle(f"disk R = {R} m, exact area πR² = {EXACT:.4f} m²", fontsize=10)
    note(fig, "Each panel has its own refinement parameter. Errors fall as the ring or the tree is refined.")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, __file__, "quadtree_vs_shoelace")


if __name__ == "__main__":
    main()
