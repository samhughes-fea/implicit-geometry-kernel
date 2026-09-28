"""
visualisation/demos/csg_expression_tree_demo.py
====================================
Geometry is a compiled program: show the CSG expression tree behind a region.

Every named region in a section is an ``Expr`` tree projected from a wall
recipe (``recipes/sdf_ir/``) and evaluated by ``sdf/ir.py``. A shear web, for
example, is literally::

    Intersect(Capsule(axis, r=t/2), Offset(Polygon(airfoil), -skin_thickness))

This demo draws that tree, evaluates it one operator at a time so the region
visibly emerges, and counts how much of the whole section's tree is *shared*
structure (the airfoil polygon is referenced by every wall but stored and
evaluated once).

Writes under ``visualisation/output/csg_expression_tree/``:

* ``csg_expression_tree.png`` - section context, node-link tree, and the
  operator-by-operator emergence of the region.
* ``csg_node_census.png`` - per-label tree size and the shared-subtree payoff.
* ``csg_expression_tree.json`` - per-label node counts, depths, unique-node total.

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/csg_expression_tree_demo.py
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
    INK,
    NAVY,
    RED,
    distance_field_norm,
    draw_distance_field,
    strip_axes_chrome,
)

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    SDFGrid,
    build_section_view,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.ir import (
    Expr,
    eval_expr_batch,
    expr_bbox,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot import (
    plot_section,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
    rotated_naca_grid_to_section_plot_axes,
)

_DPI = DEFAULT_DPI

# Nodes that combine or transform sub-fields, versus nodes that are shapes.
_OPERATOR_COLOR = "#2f6f9f"
_LEAF_COLOR = "#c2703d"
_SHARED_COLOR = "#7a4fa3"

# One line of numpy per operator, taken from the eval_expr_batch dispatch.
_OPERATOR_FORMULA = {
    "Union": "np.minimum(a, b)",
    "Intersect": "np.maximum(a, b)",
    "Subtract": "np.maximum(base, -cutter)",
    "Offset": "base - amount",
    "Rotate": "base(R(-angle) @ xy)",
    "Shell": "np.abs(base) - thickness/2",
    "Slab": "np.maximum(base, -base - depth)",
}


def _to_key(label) -> str:
    return str(getattr(label, "value", label))


def _save(fig, path: Path, *, presentation: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {path}")


def child_expressions(node: Expr) -> list[tuple[str, Expr]]:
    """Return ``(field_name, child)`` for every sub-expression of ``node``."""
    children: list[tuple[str, Expr]] = []
    for field_def in dataclasses.fields(node):
        if field_def.name.startswith("_"):
            continue
        value = getattr(node, field_def.name)
        if isinstance(value, Expr):
            children.append((field_def.name, value))
    return children


def scalar_parameters(node: Expr) -> dict[str, float]:
    """Return the non-expression scalar parameters carried by ``node``."""
    parameters: dict[str, float] = {}
    for field_def in dataclasses.fields(node):
        name = field_def.name
        if name.startswith("_") or name in ("vertices", "fn"):
            continue
        value = getattr(node, name)
        if isinstance(value, Expr):
            continue
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            parameters[name] = float(value)
    return parameters


def node_caption(node: Expr) -> str:
    """Short multi-line label for a tree node."""
    name = type(node).__name__
    if name == "Polygon":
        vertex_count = len(getattr(node, "vertices", ()))
        return f"Polygon\n{vertex_count} verts"
    if name == "CallableField":
        return f"CallableField\n{getattr(node, 'key', '')}"

    parameters = scalar_parameters(node)
    if not parameters:
        return name
    shown = list(parameters.items())[:3]
    body = "\n".join(f"{key}={value:.4g}" for key, value in shown)
    return f"{name}\n{body}"


def count_nodes(node: Expr) -> int:
    """Total node visits in the tree, counting shared subtrees each time."""
    return 1 + sum(count_nodes(child) for _, child in child_expressions(node))


def tree_depth(node: Expr) -> int:
    children = child_expressions(node)
    if not children:
        return 1
    return 1 + max(tree_depth(child) for _, child in children)


def collect_unique(node: Expr, seen: set[Expr]) -> None:
    """Accumulate structurally distinct nodes (``Expr`` is a frozen dataclass)."""
    if node in seen:
        return
    seen.add(node)
    for _, child in child_expressions(node):
        collect_unique(child, seen)


def layout_tree(root: Expr) -> tuple[dict, list]:
    """Assign (x, y) to each node: leaves in order, parents centred on children.

    Returns ``(positions, edges)`` where ``positions`` maps a node id to
    ``(x, y, node, depth)`` and ``edges`` is a list of ``(parent_id, child_id)``.
    """
    positions: dict[int, tuple[float, float, Expr, int]] = {}
    edges: list[tuple[int, int]] = []
    next_leaf_x = [0.0]

    def place(node: Expr, depth: int) -> tuple[int, float]:
        node_id = len(positions)
        children = child_expressions(node)
        # Reserve the slot before recursing so ids stay unique and stable.
        positions[node_id] = (0.0, float(-depth), node, depth)
        if not children:
            x_position = next_leaf_x[0]
            next_leaf_x[0] += 1.0
        else:
            child_x: list[float] = []
            for _, child in children:
                child_id, child_center = place(child, depth + 1)
                edges.append((node_id, child_id))
                child_x.append(child_center)
            x_position = float(np.mean(child_x))
        positions[node_id] = (x_position, float(-depth), node, depth)
        return node_id, x_position

    place(root, 0)
    return positions, edges


def draw_tree(axis, root: Expr, shared_nodes: set[Expr]) -> None:
    """Draw the expression tree as a node-link diagram.

    Leaves sit one unit apart, so a wide tree has to give up detail rather than
    overlap: captions shrink, parameter lines and the numpy annotations drop
    away once there is no room for them.
    """
    positions, edges = layout_tree(root)
    leaf_count = sum(1 for _, _, node, _ in positions.values() if not child_expressions(node))

    if leaf_count <= 6:
        box_font, show_parameters, show_formula = 10.0, True, True
    elif leaf_count <= 12:
        box_font, show_parameters, show_formula = 5.8, False, False
    else:
        box_font, show_parameters, show_formula = 4.8, False, False

    box_half_width = 0.47 if leaf_count <= 6 else 0.46
    # The numpy line rides inside the node, so the box has to be tall enough to
    # carry it. Floated outside to the right it read as stray text with nothing
    # tying it to the operator it implements.
    box_half_height = (0.42 if show_formula else 0.22) if show_parameters else 0.15
    row_pitch = 1.25 if show_parameters else 0.72

    for parent_id, child_id in edges:
        parent_x, parent_y, _, _ = positions[parent_id]
        child_x, child_y, _, _ = positions[child_id]
        axis.plot(
            [parent_x, child_x * 1.0],
            [parent_y * row_pitch, child_y * row_pitch],
            color="#9aa3ab",
            linewidth=1.3,
            zorder=1,
        )

    for node_x, node_y_raw, node, _ in positions.values():
        node_y = node_y_raw * row_pitch
        is_leaf = not child_expressions(node)
        face = _LEAF_COLOR if is_leaf else _OPERATOR_COLOR
        is_shared = node in shared_nodes
        box = FancyBboxPatch(
            (node_x - box_half_width, node_y - box_half_height),
            2.0 * box_half_width,
            2.0 * box_half_height,
            boxstyle="round,pad=0.02,rounding_size=0.06",
            facecolor=face,
            edgecolor=_SHARED_COLOR if is_shared else "white",
            linewidth=2.6 if is_shared else 1.4,
            linestyle="solid",
            zorder=2,
        )
        axis.add_patch(box)
        formula = _OPERATOR_FORMULA.get(type(node).__name__)
        carries_formula = show_formula and formula is not None
        caption = node_caption(node) if show_parameters else type(node).__name__
        axis.text(
            node_x,
            node_y + (0.08 if carries_formula else 0.0),
            caption,
            ha="center",
            va="center",
            fontsize=box_font,
            color="white",
            zorder=3,
            linespacing=1.2,
        )
        if carries_formula:
            axis.text(
                node_x,
                node_y - box_half_height + 0.09,
                formula,
                ha="center",
                va="center",
                fontsize=9.5,
                family="monospace",
                color="white",
                alpha=0.82,
                zorder=3,
            )

    all_x = [p[0] for p in positions.values()]
    all_y = [p[1] * row_pitch for p in positions.values()]
    # The formula now sits inside the node, so no right gutter is reserved for it.
    right_margin = box_half_width + 0.4
    axis.set_xlim(min(all_x) - box_half_width - 0.4, max(all_x) + right_margin)
    axis.set_ylim(
        min(all_y) - box_half_height - 0.35, max(all_y) + box_half_height + 0.35
    )
    axis.set_axis_off()


def emergence_stages(root: Expr, max_stages: int) -> list[Expr]:
    """A coherent build-up: leaf first, root last, each stage the parent of the last.

    Walks the construction spine - at each node, descend into whichever child
    carries the larger subtree, which is the region's main construction path
    rather than an incidental clip or cutter. A plain post-order traversal would
    show several unrelated subtrees and then jump to the root, which reads as a
    discontinuity rather than a build-up.
    """
    spine: list[Expr] = [root]
    node = root
    while True:
        children = [child for _, child in child_expressions(node)]
        if not children:
            break
        node = max(children, key=count_nodes)
        spine.append(node)

    spine.reverse()  # leaf -> root
    if len(spine) <= max_stages:
        return spine

    # Subsample, always keeping the leaf and the root.
    picks = np.linspace(0, len(spine) - 1, max_stages)
    indices = sorted({int(round(float(value))) for value in picks})
    indices[0] = 0
    indices[-1] = len(spine) - 1
    return [spine[index] for index in indices]


def stage_grid_for(
    root: Expr, airfoil: AirfoilSDF, nx: int, ny: int
) -> tuple[SDFGrid, tuple[float, float, float, float] | None]:
    """Grid framing the region in its local aerofoil context.

    A shear web is millimetres wide and decimetres tall, so its own bounding box
    renders as an unreadable sliver. Widen the window in x around the region so
    the aerofoil surface and the inward skin offset are both visible alongside
    it — that is what makes the ``Intersect`` legible.
    """
    bbox = expr_bbox(root)
    vertices = np.asarray(airfoil.vertices, dtype=float)
    y_min_airfoil = float(vertices[:, 1].min())
    y_max_airfoil = float(vertices[:, 1].max())
    y_pad = 0.30 * (y_max_airfoil - y_min_airfoil) + 0.01

    if bbox is None:
        x_min = float(vertices[:, 0].min()) - 0.03
        x_max = float(vertices[:, 0].max()) + 0.03
    else:
        x_center = 0.5 * (bbox[0] + bbox[1])
        # Pad by a fraction of the region width, never a multiple of it: a
        # spar cap spans most of the chord, and a 3x multiple would frame
        # several chords of empty space.
        half_width = max(
            0.16 * float(airfoil.chord), 0.62 * (bbox[1] - bbox[0])
        )
        x_min = x_center - half_width
        x_max = x_center + half_width

    grid = SDFGrid.from_bbox(
        x_min, x_max, y_min_airfoil - y_pad, y_max_airfoil + y_pad, nx=nx, ny=ny
    )
    return grid, bbox


def draw_field(
    axis,
    grid: SDFGrid,
    phi: np.ndarray,
    title: str,
    *,
    area_m2: float | None = None,
    previous_phi: np.ndarray | None = None,
    norm=None,
):
    """Signed field with the phi=0 interface picked out.

    ``previous_phi`` draws the preceding stage's interface as a dashed ghost, so
    a small ``Offset`` reads as a visible gap rather than an identical picture.

    Returns the mesh so the caller can hang one shared colourbar off the row:
    without it the three panels show magnitude with no scale to read it against.
    """
    mesh = draw_distance_field(axis, grid.X, grid.Y, phi, norm=norm, zero_line_width=1.5)
    if previous_phi is not None:
        axis.contour(
            grid.X,
            grid.Y,
            previous_phi,
            levels=[0.0],
            colors="#5b6672",
            linewidths=1.1,
            linestyles="dashed",
        )
    axis.set_title(title, fontsize=13, loc="left", color=INK)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_aspect("equal")
    if area_m2 is not None:
        axis.text(
            0.03,
            0.04,
            f"enclosed {area_m2 * 1e4:.1f} cm²",
            transform=axis.transAxes,
            fontsize=11,
            color="#2b3238",
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 1.8},
        )
    return mesh


def _annotate_target(axis, x_plot, y_plot, target_phi, label: str) -> None:
    """Leader-line the highlighted region so it is named where it is drawn.

    The colour alone identifies it only for someone who already knows which
    region is which; from the back of a room the label has to be on the figure.
    """
    inside = np.asarray(target_phi) < 0.0
    if not inside.any():
        return
    x_center = float(np.asarray(x_plot)[inside].mean())
    y_center = float(np.asarray(y_plot)[inside].mean())

    _, y_top = axis.get_ylim()
    axis.annotate(
        label,
        xy=(x_center, y_center),
        xytext=(x_center, y_top * 0.72),
        ha="center",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        color=INK,
        zorder=6,
        arrowprops={"arrowstyle": "-", "color": INK, "linewidth": 0.9},
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--key", default="2D-F")
    parser.add_argument("--naca", default="2412")
    parser.add_argument("--chord", type=float, default=1.0)
    parser.add_argument("--twist-deg", type=float, default=0.0)
    parser.add_argument("--nx", type=int, default=900)
    parser.add_argument("--ny", type=int, default=400)
    parser.add_argument(
        "--label",
        default="shear_web_0",
        help="Component whose expression tree is drawn.",
    )
    parser.add_argument("--skin-thickness", type=float, default=0.0035)
    parser.add_argument("--spar-cap-thickness", type=float, default=0.012)
    parser.add_argument(
        "--shear-web-thickness",
        type=float,
        default=0.008,
        help="Slightly generous by default so the web reads at slide size.",
    )
    parser.add_argument("--max-stages", type=int, default=4)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--presentation", action="store_true", help="Deck styling.")
    args = parser.parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)

    out_dir = args.out_dir or (
        Path(__file__).resolve().parents[1] / "output" / "csg_expression_tree"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = resolve_system_type(str(args.key))
    twist_rad = float(np.deg2rad(float(args.twist_deg)))
    airfoil = AirfoilSDF.from_naca(str(args.naca), chord=float(args.chord), n_points=180)
    section = build_section_view(
        airfoil,
        layout,
        twist_rad=twist_rad,
        skin_thickness=float(args.skin_thickness),
        spar_cap_thickness=float(args.spar_cap_thickness),
        shear_web_thickness=float(args.shear_web_thickness),
    )

    label_keys = [_to_key(label) for label in section]
    if args.label not in label_keys:
        raise SystemExit(
            f"--label {args.label!r} not in section {args.key}; available: {label_keys}"
        )

    airfoil_for_grid = airfoil.rotate(twist_rad) if abs(twist_rad) > 1e-12 else airfoil
    section_grid = SDFGrid.from_airfoil(
        airfoil_for_grid, nx=int(args.nx), ny=int(args.ny), padding=0.08
    )

    target_expression = section[args.label].expr

    # --- census over the whole section -------------------------------------
    per_label: dict[str, dict[str, int]] = {}
    all_unique: set[Expr] = set()
    total_visits = 0
    for label in section:
        key = _to_key(label)
        expression = getattr(section[label], "expr", None)
        if expression is None:
            continue
        visits = count_nodes(expression)
        own_unique: set[Expr] = set()
        collect_unique(expression, own_unique)
        collect_unique(expression, all_unique)
        total_visits += visits
        per_label[key] = {
            "node_visits": int(visits),
            "unique_nodes": int(len(own_unique)),
            "depth": int(tree_depth(expression)),
        }

    shared_nodes: set[Expr] = set()
    for candidate in all_unique:
        appearances = 0
        for label in section:
            expression = getattr(section[label], "expr", None)
            if expression is None:
                continue
            own: set[Expr] = set()
            collect_unique(expression, own)
            if candidate in own:
                appearances += 1
                if appearances > 1:
                    break
        if appearances > 1:
            shared_nodes.add(candidate)

    # --- figure 1: tree + emergence ----------------------------------------
    stages = emergence_stages(target_expression, int(args.max_stages))
    stage_grid, target_bbox = stage_grid_for(
        target_expression, airfoil_for_grid, nx=460, ny=300
    )

    figure = plt.figure(figsize=(13.5, 7.6))
    # A station is a thin lens. Given a tall top row, equal aspect fills the
    # difference with dead space above and below the aerofoil rather than with
    # geometry, so the row is deliberately short.
    outer = figure.add_gridspec(2, 1, height_ratios=[0.82, 1.0], hspace=0.22)
    # The tree is three or four nodes and does not need the larger half; the
    # section is a thin lens and every millimetre of width becomes drawn chord.
    top = outer[0].subgridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.10)
    axis_section = figure.add_subplot(top[0, 0])
    axis_tree = figure.add_subplot(top[0, 1])

    # The assembly is context, not subject: hold it at low alpha so the one
    # region this figure is about can be filled solid over the top. The legend
    # goes with it — the anatomy figure names all nine regions one slide
    # earlier, and repeating them here costs the section its best real estate.
    plot_section(
        section,
        section_grid,
        ax=axis_section,
        alpha=0.16,
        show_airfoil=False,
        show_legend=False,
        auto_limits=True,
        title="",
    )
    strip_axes_chrome(axis_section)
    axis_section.set_title(f"{args.key}: the region built below is {args.label}", loc="left",
                           fontsize=13, color=INK)
    # Pick the demonstrated region out of the assembly. plot_section draws in
    # the structural frame (horizontal = Y_B, vertical = X_B), so the SDF-grid
    # arrays are passed through the same remap before contouring.
    target_phi_section = np.asarray(
        eval_expr_batch(target_expression, section_grid.X, section_grid.Y, {}),
        dtype=float,
    )
    x_plot, y_plot = rotated_naca_grid_to_section_plot_axes(
        section_grid.X, section_grid.Y, twist_rad=twist_rad, chord_m=float(args.chord)
    )
    # Filled, not merely outlined. A 3.5 mm web drawn as a 2.4 pt contour at
    # section scale is indistinguishable from a plotting artefact — the earlier
    # version promised a highlight and delivered a thin black rectangle.
    axis_section.contourf(
        x_plot,
        y_plot,
        target_phi_section,
        levels=[float(np.nanmin(target_phi_section)) - 1.0, 0.0],
        colors=[RED],
        alpha=0.95,
        zorder=4,
    )
    axis_section.contour(
        x_plot,
        y_plot,
        target_phi_section,
        levels=[0.0],
        colors=[NAVY],
        linewidths=2.0,
        zorder=5,
    )
    _annotate_target(axis_section, x_plot, y_plot, target_phi_section, args.label)

    axis_tree.set_title(
        "Its expression tree · blue operators, orange leaves,\n"
        "purple outline = shared with another region",
        loc="left", fontsize=12, color=INK,
    )
    draw_tree(axis_tree, target_expression, shared_nodes)

    bottom = outer[1].subgridspec(1, len(stages), wspace=0.07)
    stage_records: list[dict] = []
    previous_phi: np.ndarray | None = None
    stage_axes: list = []
    stage_mesh = None

    # A fresh memo per stage, deliberately. The memo is keyed by Expr alone, but
    # a Rotate node evaluates its base at pulled-back coordinates - so a memo
    # carried across a Rotate boundary would hand back values sampled at the
    # wrong points. Sharing is safe within one root (as the package does); it is
    # not safe across these independently rooted stages.
    stage_fields = [
        np.asarray(
            eval_expr_batch(stage_expression, stage_grid.X, stage_grid.Y, {}),
            dtype=float,
        )
        for stage_expression in stages
    ]

    # One norm for the row, taken from the first stage. Per-panel norms make the
    # progression unreadable: the final Intersect encloses a 3.5 mm web, so a
    # norm derived from it clips everything past 3.5 mm and floods the panel
    # with one flat colour. The stages are the same field being cut down — they
    # have to be measured against the same scale to show that.
    stage_norm = distance_field_norm(stage_fields[0])

    for stage_index, stage_expression in enumerate(stages):
        axis_stage = figure.add_subplot(bottom[0, stage_index])
        phi = stage_fields[stage_index]
        # Ghost the previous interface only where this stage transformed it
        # in place (Offset / Shell / Slab); a fresh leaf has nothing to compare.
        ghost = (
            previous_phi
            if type(stage_expression).__name__ in ("Offset", "Shell", "Slab")
            else None
        )
        stage_area = stage_grid.area(phi)
        stage_mesh = draw_field(
            axis_stage,
            stage_grid,
            phi,
            f"{stage_index + 1}. {type(stage_expression).__name__}",
            area_m2=stage_area,
            previous_phi=ghost,
            norm=stage_norm,
        )
        stage_axes.append(axis_stage)
        stage_records.append(
            {
                "node": type(stage_expression).__name__,
                "enclosed_area_m2": float(stage_area),
            }
        )
        previous_phi = phi

    # One bar for the row: the stage panels carry magnitude and had no scale to
    # read it against. All three share ``stage_norm``, so one bar is honest.
    if stage_mesh is not None and stage_axes:
        colour_bar = figure.colorbar(
            stage_mesh, ax=stage_axes, orientation="vertical", fraction=0.022, pad=0.012
        )
        colour_bar.set_label("φ [m]", fontsize=11)
        colour_bar.ax.tick_params(labelsize=10)
        colour_bar.outline.set_visible(False)

    _save(figure, out_dir / "csg_expression_tree.png", presentation=presentation)

    # --- figure 2: node census ---------------------------------------------
    census_figure, census_axes = plt.subplots(1, 2, figsize=(15.0, 5.6))
    axis_bars, axis_summary = census_axes

    ordered_labels = sorted(per_label, key=lambda k: per_label[k]["node_visits"], reverse=True)
    y_positions = np.arange(len(ordered_labels))
    visit_values = [per_label[k]["node_visits"] for k in ordered_labels]
    unique_values = [per_label[k]["unique_nodes"] for k in ordered_labels]

    axis_bars.barh(
        y_positions, visit_values, height=0.62, color=_OPERATOR_COLOR, label="node visits"
    )
    axis_bars.barh(
        y_positions, unique_values, height=0.32, color=_LEAF_COLOR, label="distinct nodes"
    )
    axis_bars.set_yticks(y_positions)
    axis_bars.set_yticklabels(ordered_labels, fontsize=8.5)
    axis_bars.invert_yaxis()
    axis_bars.set_xlabel("nodes in the region's expression tree")
    axis_bars.set_title("Per-region tree size", fontsize=11)
    axis_bars.legend(fontsize=8.5)
    axis_bars.grid(axis="x", alpha=0.25)

    unique_total = len(all_unique)
    sharing_factor = (total_visits / unique_total) if unique_total else float("nan")
    axis_summary.bar(
        ["node visits\n(naive walk)", "distinct nodes\n(memoised)"],
        [total_visits, unique_total],
        color=[_OPERATOR_COLOR, _SHARED_COLOR],
        width=0.55,
    )
    for index, value in enumerate((total_visits, unique_total)):
        axis_summary.text(
            index, value, f"{value}", ha="center", va="bottom", fontsize=12, fontweight="bold"
        )
    axis_summary.set_title(
        f"Whole section {args.key}: shared-subtree memoisation\n"
        f"{sharing_factor:.2f}x fewer field evaluations",
        fontsize=11,
    )
    axis_summary.set_ylabel("nodes")
    axis_summary.grid(axis="y", alpha=0.25)
    axis_summary.text(
        0.74,
        0.60,
        f"{len(shared_nodes)} of the {unique_total} distinct nodes\n"
        "are referenced by more than one\nregion - the aerofoil polygon\namong them",
        transform=axis_summary.transAxes,
        ha="center",
        va="center",
        fontsize=9.5,
        color="#41505c",
        bbox={
            "facecolor": "white",
            "alpha": 0.9,
            "edgecolor": "#d3d8dd",
            "boxstyle": "round,pad=0.5",
        },
    )

    census_figure.suptitle(
        "The section's regions are distinct trees over a shared pool of nodes",
        fontsize=13,
    )
    census_figure.tight_layout()
    _save(census_figure, out_dir / "csg_node_census.png")

    payload = {
        "key": str(args.key),
        "naca": str(args.naca),
        "chord_m": float(args.chord),
        "twist_deg": float(args.twist_deg),
        "label": str(args.label),
        "thicknesses_m": {
            "skin": float(args.skin_thickness),
            "spar_cap": float(args.spar_cap_thickness),
            "shear_web": float(args.shear_web_thickness),
        },
        "per_label": per_label,
        "total_node_visits": int(total_visits),
        "distinct_nodes": int(unique_total),
        "sharing_factor": float(sharing_factor),
        "nodes_shared_across_regions": int(len(shared_nodes)),
        "demonstrated_tree": {
            "root": type(target_expression).__name__,
            "depth": int(tree_depth(target_expression)),
            "node_visits": int(count_nodes(target_expression)),
            "stages": stage_records,
        },
    }
    json_path = out_dir / "csg_expression_tree.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
