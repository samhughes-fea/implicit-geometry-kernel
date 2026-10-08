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
* ``csg_expression_tree_stages.gif`` - looping dwell on each emergence stage
  (Polygon → Offset → Intersect for a shear web). Deck live media; the PNG is
  the PDF still.
* ``csg_node_census.png`` - per-label tree size and the shared-subtree payoff.
* ``csg_expression_tree.json`` - per-label node counts, depths, unique-node total.

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/csg_expression_tree_demo.py --presentation
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
    LIGHT_TEAL,
    MUTED,
    NAVY,
    RED,
    SURFACE,
    TEAL,
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
    Capsule,
    Expr,
    Offset,
    Polygon,
    Slab,
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

# Deck-matched node fill: navy operators (hdr), surface leaves (chain-box),
# teal type labels + rim for shared — same tokens as the HTML deck :root.
_OPERATOR_FACE = NAVY
_OPERATOR_EDGE = NAVY
_LEAF_FACE = SURFACE
_LEAF_EDGE = TEAL
_SHARED_EDGE = TEAL
_CONNECTOR = MUTED
_OPERATOR_TEXT = "#ffffff"
_LEAF_TEXT = MUTED
_FORMULA_OPERATOR_ALPHA = 0.78

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


def _write_stages_gif(
    out_path: Path,
    *,
    stages: list[Expr],
    stage_fields: list[np.ndarray],
    stage_grid: SDFGrid,
    stage_norm,
    label: str,
    dwell_ms: int = 1100,
) -> None:
    """Looping GIF: one frame per emergence stage, same field scale as the PNG row."""
    from io import BytesIO

    from PIL import Image

    frames: list[Image.Image] = []
    previous_phi: np.ndarray | None = None
    names = [type(stage).__name__ for stage in stages]

    for stage_index, (stage_expression, phi) in enumerate(zip(stages, stage_fields)):
        figure, axis = plt.subplots(figsize=(7.2, 5.2))
        ghost = (
            previous_phi
            if type(stage_expression).__name__ in ("Offset", "Shell", "Slab")
            else None
        )
        draw_field(
            axis,
            stage_grid,
            phi,
            f"{stage_index + 1}. {type(stage_expression).__name__}",
            area_m2=stage_grid.area(phi),
            previous_phi=ghost,
            norm=stage_norm,
        )
        progress = "  →  ".join(
            f"[{name}]" if index == stage_index else name
            for index, name in enumerate(names)
        )
        figure.suptitle(
            f"{label} · {progress}",
            fontsize=12,
            color=INK,
            x=0.02,
            ha="left",
            y=0.98,
        )
        figure.subplots_adjust(left=0.04, right=0.98, top=0.90, bottom=0.04)
        buffer = BytesIO()
        figure.savefig(buffer, format="png", dpi=120, facecolor="white")
        plt.close(figure)
        buffer.seek(0)
        frames.append(Image.open(buffer).convert("P", palette=Image.ADAPTIVE))
        previous_phi = phi

    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        out_path,
        save_all=True,
        append_images=frames[1:],
        duration=dwell_ms,
        loop=0,
        optimize=True,
    )
    print(f"Wrote {out_path}")


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
    if name == "Capsule":
        # Two lines only — long endpoint dumps clip inside the box.
        parameters = scalar_parameters(node)
        radius = parameters.get("r")
        ax_ = parameters.get("ax")
        bx_ = parameters.get("bx")
        if radius is not None and ax_ is not None and bx_ is not None:
            if abs(ax_ - bx_) < 1e-9:
                return f"Capsule\nx={ax_:.3g} · r={radius:.3g}"
            return f"Capsule\nr={radius:.3g}"
        if radius is not None:
            return f"Capsule\nr={radius:.3g}"
        return "Capsule"

    parameters = scalar_parameters(node)
    if not parameters:
        return name
    if name == "Offset" and "amount" in parameters:
        return f"Offset\namount={parameters['amount']:.4g}"
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


def layout_tree(root: Expr, *, leaf_gap: float = 1.0) -> tuple[dict, list]:
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
            next_leaf_x[0] += leaf_gap
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


def draw_tree(
    axis,
    root: Expr,
    shared_nodes: set[Expr],
    *,
    names_only: bool = False,
) -> None:
    """Draw the expression tree in deck chrome (navy / teal / surface).

    Matches the CDFAM HTML chain-box / gap-card language rather than a generic
    orange–blue diagram palette. Small trees get wide leaf gaps and elbows;
    wide trees shrink type and drop parameter / formula lines rather than overlap.
    ``names_only`` keeps type names alone (deck wall panels).
    """
    leaf_count = sum(1 for _ in _iter_leaves(root))
    leaf_gap = 2.15 if leaf_count <= 3 else (1.55 if leaf_count <= 6 else 1.0)
    positions, edges = layout_tree(root, leaf_gap=leaf_gap)

    if names_only:
        box_font, formula_font, show_parameters, show_formula = 13.0, 9.0, False, False
        box_half_width, box_half_height, row_pitch = 0.85, 0.28, 1.15
    elif leaf_count <= 3:
        box_font, formula_font, show_parameters, show_formula = 11.5, 9.0, True, True
        box_half_width, box_half_height, row_pitch = 0.92, 0.50, 1.65
    elif leaf_count <= 6:
        box_font, formula_font, show_parameters, show_formula = 10.0, 8.0, True, True
        box_half_width, box_half_height, row_pitch = 0.70, 0.44, 1.40
    elif leaf_count <= 12:
        box_font, formula_font, show_parameters, show_formula = 5.8, 5.0, False, False
        box_half_width, box_half_height, row_pitch = 0.46, 0.15, 0.72
    else:
        box_font, formula_font, show_parameters, show_formula = 4.8, 4.5, False, False
        box_half_width, box_half_height, row_pitch = 0.46, 0.15, 0.72

    for parent_id, child_id in edges:
        parent_x, parent_y_raw, _, _ = positions[parent_id]
        child_x, child_y_raw, _, _ = positions[child_id]
        parent_y = parent_y_raw * row_pitch
        child_y = child_y_raw * row_pitch
        mid_y = 0.5 * (parent_y + child_y)
        axis.plot(
            [parent_x, parent_x, child_x, child_x],
            [parent_y - box_half_height, mid_y, mid_y, child_y + box_half_height],
            color=_CONNECTOR,
            linewidth=1.5,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=1,
        )

    for node_x, node_y_raw, node, _ in positions.values():
        node_y = node_y_raw * row_pitch
        is_leaf = not child_expressions(node)
        is_shared = node in shared_nodes
        face = _LEAF_FACE if is_leaf else _OPERATOR_FACE
        if is_shared:
            edge, edge_width = _SHARED_EDGE, 2.4
        elif is_leaf:
            edge, edge_width = _LEAF_EDGE, 1.6
        else:
            edge, edge_width = _OPERATOR_EDGE, 1.2
        box = FancyBboxPatch(
            (node_x - box_half_width, node_y - box_half_height),
            2.0 * box_half_width,
            2.0 * box_half_height,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            facecolor=face,
            edgecolor=edge,
            linewidth=edge_width,
            linestyle="solid",
            zorder=2,
        )
        axis.add_patch(box)
        formula = _OPERATOR_FORMULA.get(type(node).__name__)
        carries_formula = show_formula and formula is not None and not is_leaf
        caption = node_caption(node) if show_parameters else type(node).__name__
        # First line = type name (teal weight on deck cards); rest = params.
        caption_lines = caption.split("\n")
        name_line = caption_lines[0]
        param_lines = caption_lines[1:]
        text_color = _LEAF_TEXT if is_leaf else _OPERATOR_TEXT
        name_y = node_y + (0.12 if (param_lines or carries_formula) else 0.0)
        axis.text(
            node_x,
            name_y,
            name_line,
            ha="center",
            va="center",
            fontsize=box_font,
            color=TEAL if is_leaf else _OPERATOR_TEXT,
            zorder=3,
            fontweight="bold",
        )
        if param_lines:
            axis.text(
                node_x,
                node_y - (0.06 if carries_formula else 0.10),
                "\n".join(param_lines),
                ha="center",
                va="center",
                fontsize=box_font - 1.5,
                color=text_color,
                zorder=3,
                linespacing=1.2,
            )
        if carries_formula:
            axis.text(
                node_x,
                node_y - box_half_height + 0.11,
                formula,
                ha="center",
                va="center",
                fontsize=formula_font,
                family="monospace",
                color=_OPERATOR_TEXT,
                alpha=_FORMULA_OPERATOR_ALPHA,
                zorder=3,
            )

    all_x = [p[0] for p in positions.values()]
    all_y = [p[1] * row_pitch for p in positions.values()]
    right_margin = box_half_width + 0.35
    axis.set_xlim(min(all_x) - box_half_width - 0.35, max(all_x) + right_margin)
    axis.set_ylim(
        min(all_y) - box_half_height - 0.30, max(all_y) + box_half_height + 0.30
    )
    axis.set_facecolor(SURFACE)
    axis.set_axis_off()


def _iter_leaves(node: Expr):
    children = child_expressions(node)
    if not children:
        yield node
        return
    for _, child in children:
        yield from _iter_leaves(child)


def _airfoil_polygon(airfoil) -> Polygon:
    verts = np.asarray(airfoil.vertices, dtype=float)
    return Polygon(tuple((float(x), float(y)) for x, y in verts))


def pedagogical_wall_expr(
    wall: str,
    *,
    airfoil,
    skin_thickness: float,
    spar_cap_thickness: float,
    shear_web_thickness: float,
    production_web: Expr | None = None,
) -> Expr:
    """Compact recipe spine for deck wall panels (not the full winner-map tree).

    Skin/cap production exprs carry half-planes and web subtracts that overwhelm
    a slide; the talk track is Offset / Slab / Capsule∩cavity.
    """
    leaf = _airfoil_polygon(airfoil)
    cavity = Offset(leaf, -float(skin_thickness))
    kind = str(wall).strip().lower()
    if kind in {"skin", "skin_suction_surface"}:
        return cavity
    if kind in {"cap", "spar_cap", "spar_cap_suction_surface"}:
        return Slab(cavity, float(spar_cap_thickness))
    if kind in {"web", "shear_web", "shear_web_0"}:
        if production_web is not None:
            return production_web
        # Fallback axis near mid-chord if the section build was skipped.
        half = 0.5 * float(shear_web_thickness)
        return Capsule(0.35, -0.08, 0.35, 0.08, half)
    raise ValueError(f"unknown wall {wall!r}; expected skin|cap|web")


def write_compact_tree_png(
    path: Path,
    root: Expr,
    *,
    shared_nodes: set[Expr] | None = None,
    title: str,
    presentation: bool = False,
) -> None:
    """Tree-only panel for the wall trio (no section context, no emergence row)."""
    shared = shared_nodes if shared_nodes is not None else set()
    figure, axis = plt.subplots(figsize=(4.6, 4.2))
    figure.patch.set_facecolor("white")
    axis.set_title(title, loc="left", fontsize=11, color=MUTED, pad=8)
    draw_tree(axis, root, shared, names_only=True)
    figure.subplots_adjust(left=0.04, right=0.96, top=0.88, bottom=0.04)
    _save(figure, path, presentation=presentation)

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
    parser.add_argument(
        "--names-only",
        action="store_true",
        help="Type names on tree nodes (no parameters / formulas).",
    )
    parser.add_argument(
        "--tree-only",
        action="store_true",
        help="Write compact tree-only PNG(s); skip section/emergence/census figures.",
    )
    parser.add_argument(
        "--deck-walls",
        action="store_true",
        help="Write csg_tree_skin/spar_cap/shear_web name-only panels for slides 13–15.",
    )
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
    if args.label not in label_keys and not args.deck_walls:
        raise SystemExit(
            f"--label {args.label!r} not in section {args.key}; available: {label_keys}"
        )

    if args.deck_walls:
        web_expr = section["shear_web_0"].expr
        panels = (
            (
                "skin",
                out_dir / "csg_tree_skin.png",
                "Skin · Offset of the NACA leaf",
            ),
            (
                "cap",
                out_dir / "csg_tree_spar_cap.png",
                "Spar cap · Slab on the cavity",
            ),
            (
                "web",
                out_dir / "csg_tree_shear_web.png",
                "Shear web · Capsule ∩ cavity",
            ),
        )
        for wall, path, title in panels:
            expr = pedagogical_wall_expr(
                wall,
                airfoil=airfoil,
                skin_thickness=float(args.skin_thickness),
                spar_cap_thickness=float(args.spar_cap_thickness),
                shear_web_thickness=float(args.shear_web_thickness),
                production_web=web_expr if wall == "web" else None,
            )
            write_compact_tree_png(
                path,
                expr,
                shared_nodes=set(),
                title=title,
                presentation=presentation,
            )
        if args.tree_only:
            return 0


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
    # Tree gets equal width to the section so nodes and elbow connectors breathe.
    top = outer[0].subgridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.12)
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
        "Expression tree · navy = operators · teal = leaves\n"
        "teal rim = shared subtree",
        loc="left", fontsize=12, color=MUTED,
    )
    draw_tree(
        axis_tree,
        target_expression,
        shared_nodes,
        names_only=bool(args.names_only),
    )

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

    _write_stages_gif(
        out_dir / "csg_expression_tree_stages.gif",
        stages=stages,
        stage_fields=stage_fields,
        stage_grid=stage_grid,
        stage_norm=stage_norm,
        label=str(args.label),
    )

    # --- figure 2: node census ---------------------------------------------
    census_figure, census_axes = plt.subplots(1, 2, figsize=(15.0, 5.6))
    axis_bars, axis_summary = census_axes

    ordered_labels = sorted(per_label, key=lambda k: per_label[k]["node_visits"], reverse=True)
    y_positions = np.arange(len(ordered_labels))
    visit_values = [per_label[k]["node_visits"] for k in ordered_labels]
    unique_values = [per_label[k]["unique_nodes"] for k in ordered_labels]

    axis_bars.barh(
        y_positions, visit_values, height=0.62, color=NAVY, label="node visits"
    )
    axis_bars.barh(
        y_positions, unique_values, height=0.32, color=TEAL, label="distinct nodes"
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
        color=[NAVY, TEAL],
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
