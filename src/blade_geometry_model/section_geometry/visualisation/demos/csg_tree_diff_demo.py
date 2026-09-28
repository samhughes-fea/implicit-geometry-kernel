"""
visualisation/demos/csg_tree_diff_demo.py
=========================================
Diff two design vectors by comparing their CSG expression trees.

Builds two sections (same layout, one gauge moved by default) and draws one
merged tree for a chosen region: each node is reused as-is, rebuilt because a
child changed, or re-parameterised (shown as ``A -> B``). Both sections are
drawn with the affected regions highlighted, and the whole-section census of
distinct nodes goes in the legend. Reuses layout_tree helpers from
csg_expression_tree_demo.py.

Writes under visualisation/output/csg_tree_diff/:

* csg_tree_diff.png
* csg_tree_diff.json

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/csg_tree_diff_demo.py --presentation
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

_tree = importlib.util.spec_from_file_location(
    "csg_expression_tree_demo",
    Path(__file__).resolve().parent / "csg_expression_tree_demo.py",
)
_tree_mod = importlib.util.module_from_spec(_tree)
assert _tree.loader is not None
_tree.loader.exec_module(_tree_mod)

child_expressions = _tree_mod.child_expressions
collect_unique = _tree_mod.collect_unique
count_nodes = _tree_mod.count_nodes
layout_tree = _tree_mod.layout_tree
node_caption = _tree_mod.node_caption
scalar_parameters = _tree_mod.scalar_parameters

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    build_section_view,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.ir import Expr
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    add_output_and_presentation_args,
    apply_presentation_style,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
    INK,
    MUTED,
    NAVY,
    RED,
    strip_axes_chrome,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot import (
    plot_section_hybrid,
)

_REUSED = "#9aa3ae"
_REBUILT = "#f3b7c4"
_CHANGED = RED
_UNTOUCHED_REGION = "#d9dde3"
_AFFECTED_REGION = "#e38aa0"


def _to_key(label) -> str:
    return str(getattr(label, "value", label))


def _node_status(node_a: Expr, node_b: Expr) -> str:
    """How node B relates to node A at the same place in the tree.

    ``reused``: the whole subtree is identical. ``changed``: this node's own
    parameters moved. ``rebuilt``: its own parameters are equal, but something
    beneath it changed.
    """
    # Expr is a frozen dataclass: equality is structural, vertices included —
    # the same notion of "distinct node" collect_unique counts on slide 9.
    if node_a == node_b:
        return "reused"
    if type(node_a) is not type(node_b):
        raise SystemExit(
            f"Trees differ in structure ({type(node_a).__name__} vs {type(node_b).__name__}); "
            "this demo diffs parameter changes on one layout."
        )
    if scalar_parameters(node_a) != scalar_parameters(node_b):
        return "changed"
    return "rebuilt"


def _paired_walk(node_a: Expr, node_b: Expr):
    """Yield ``(node_a, node_b)`` in the same pre-order ``layout_tree`` numbers nodes."""
    yield node_a, node_b
    children_a = child_expressions(node_a)
    children_b = child_expressions(node_b)
    if len(children_a) != len(children_b):
        raise SystemExit("Trees differ in arity; this demo diffs parameter changes on one layout.")
    for (_, child_a), (_, child_b) in zip(children_a, children_b):
        yield from _paired_walk(child_a, child_b)


def _format_length(value: float) -> str:
    return f"{value * 1e3:.1f} mm" if abs(value) < 1.0 else f"{value:.3g}"


def _changed_parameter_lines(node_a: Expr, node_b: Expr, limit: int = 2) -> list[str]:
    parameters_a = scalar_parameters(node_a)
    parameters_b = scalar_parameters(node_b)
    lines = [
        f"{name}  {_format_length(parameters_a[name])} → {_format_length(parameters_b[name])}"
        for name in parameters_a
        if name in parameters_b and parameters_a[name] != parameters_b[name]
    ]
    return lines[:limit]


def _draw_merged_tree(axis, root_a: Expr, root_b: Expr) -> None:
    positions, edges = layout_tree(root_a)
    pairs = list(_paired_walk(root_a, root_b))
    for parent_id, child_id in edges:
        x0, y0, _, _ = positions[parent_id]
        x1, y1, _, _ = positions[child_id]
        axis.plot([x0, x1], [y0, y1], color=MUTED, lw=1.2, zorder=1)

    fills = {"reused": _REUSED, "rebuilt": _REBUILT, "changed": _CHANGED}
    text_colours = {"reused": "white", "rebuilt": INK, "changed": "white"}
    details = {"reused": ["reused as-is"], "rebuilt": ["rebuilt: a child changed"]}
    for node_id, (x_position, y_position, _, _) in positions.items():
        node_a, node_b = pairs[node_id]
        status = _node_status(node_a, node_b)
        detail = details.get(status) or _changed_parameter_lines(node_a, node_b)
        box = FancyBboxPatch(
            (x_position - 0.46, y_position - 0.24),
            0.92,
            0.48,
            boxstyle="round,pad=0.02,rounding_size=0.06",
            linewidth=0.0,
            facecolor=fills[status],
            zorder=2,
        )
        axis.add_patch(box)
        axis.text(
            x_position, y_position + 0.10, type(node_a).__name__, ha="center", va="center",
            fontsize=13, fontweight="bold", color=text_colours[status], zorder=3,
        )
        axis.text(
            x_position, y_position - 0.10, "\n".join(detail), ha="center", va="center",
            fontsize=10, color=text_colours[status], zorder=3, linespacing=1.15,
        )
    axis.set_axis_off()
    xs = [position[0] for position in positions.values()]
    ys = [position[1] for position in positions.values()]
    axis.set_xlim(min(xs) - 0.6, max(xs) + 0.6)
    axis.set_ylim(min(ys) - 0.45, max(ys) + 0.45)


def _region_expr(section, label: str) -> Expr:
    label_keys = [_to_key(key) for key in section]
    if label not in label_keys:
        raise SystemExit(f"Label {label!r} not found; have {label_keys}")
    return section[label].expr


def _section_diff(section_a, section_b) -> dict:
    """Whole-section census: every distinct node of B, by how it relates to A."""
    statuses: dict[Expr, str] = {}
    changed_regions: set[str] = set()
    for label in section_a:
        key = _to_key(label)
        root_a = _region_expr(section_a, key)
        root_b = _region_expr(section_b, key)
        for node_a, node_b in _paired_walk(root_a, root_b):
            statuses.setdefault(node_b, _node_status(node_a, node_b))
        if _node_status(root_a, root_b) != "reused":
            changed_regions.add(key)
    counts = {
        status: sum(1 for value in statuses.values() if value == status)
        for status in ("reused", "rebuilt", "changed")
    }
    return {
        "counts": counts,
        "changed_regions": sorted(changed_regions),
        "region_count": len(list(section_a)),
    }


def _draw_section(axis, section, changed_regions: set[str], label: str, title: str) -> None:
    def colour(region_label):
        key = _to_key(region_label)
        if key == label:
            return RED
        return _AFFECTED_REGION if key in changed_regions else _UNTOUCHED_REGION

    plot_section_hybrid(
        section,
        residual_mesh=None,
        ax=axis,
        alpha=0.9,
        core_alpha=0.35,
        title="",
        show_legend=False,
        color_fn=colour,
    )
    axis.set_aspect("equal", adjustable="datalim")
    strip_axes_chrome(axis)
    axis.set_title(title, fontsize=13, color=INK, loc="left")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="csg_tree_diff")
    parser.add_argument("--layout", default="2D-F")
    parser.add_argument("--label", default="shear_web_0")
    parser.add_argument("--skin-a", type=float, default=0.004)
    parser.add_argument("--skin-b", type=float, default=0.004)
    parser.add_argument("--cap-a", type=float, default=0.012)
    parser.add_argument("--cap-b", type=float, default=0.012)
    parser.add_argument("--web-a", type=float, default=0.005)
    parser.add_argument("--web-b", type=float, default=0.008)
    args = parser.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    airfoil = AirfoilSDF.from_naca("2412", chord=1.0, n_points=160)
    layout = resolve_system_type(str(args.layout))

    def build(skin: float, cap: float, web: float):
        return build_section_view(
            airfoil,
            layout,
            skin_thickness=float(skin),
            spar_cap_thickness=float(cap),
            shear_web_thickness=float(web),
        )

    section_a = build(args.skin_a, args.cap_a, args.web_a)
    section_b = build(args.skin_b, args.cap_b, args.web_b)
    left = _region_expr(section_a, str(args.label))
    right = _region_expr(section_b, str(args.label))
    census = _section_diff(section_a, section_b)
    changed_regions = set(census["changed_regions"])

    moved = [
        f"{name} {value_a * 1e3:g} → {value_b * 1e3:g} mm"
        for name, value_a, value_b in (
            ("skin", args.skin_a, args.skin_b),
            ("cap", args.cap_a, args.cap_b),
            ("web", args.web_a, args.web_b),
        )
        if value_a != value_b
    ]

    fig = plt.figure(figsize=(14.0, 5.8))
    grid = fig.add_gridspec(2, 2, width_ratios=(1.2, 1.0), wspace=0.04, hspace=0.25)
    _draw_section(fig.add_subplot(grid[0, 0]), section_a, set(), "", f"A · {args.layout}")
    _draw_section(
        fig.add_subplot(grid[1, 0]),
        section_b,
        changed_regions,
        str(args.label),
        f"B · {', '.join(moved) or 'no change'}",
    )
    tree_axis = fig.add_subplot(grid[:, 1])
    _draw_merged_tree(tree_axis, left, right)
    tree_axis.set_title(f"{args.label}: one tree, A → B", fontsize=13, color=INK)

    counts = census["counts"]
    total = sum(counts.values())
    handles = [
        FancyBboxPatch((0, 0), 1, 1, facecolor=_CHANGED, linewidth=0,
                       label=f"parameter changed ({counts['changed']})"),
        FancyBboxPatch((0, 0), 1, 1, facecolor=_REBUILT, linewidth=0,
                       label=f"rebuilt above a change ({counts['rebuilt']})"),
        FancyBboxPatch((0, 0), 1, 1, facecolor=_REUSED, linewidth=0,
                       label=f"reused unchanged ({counts['reused']})"),
    ]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.0),
        ncol=3,
        frameon=False,
        fontsize=12,
        title=f"whole section: {total} distinct nodes in B",
        title_fontsize=12,
        alignment="left",
    )

    png_path = out_dir / "csg_tree_diff.png"
    save_figure(fig, png_path, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {png_path}")

    left_unique: set[Expr] = set()
    right_unique: set[Expr] = set()
    collect_unique(left, left_unique)
    collect_unique(right, right_unique)
    payload = {
        "layout": str(args.layout),
        "label": str(args.label),
        "design_a": {
            "skin_thickness": float(args.skin_a),
            "spar_cap_thickness": float(args.cap_a),
            "shear_web_thickness": float(args.web_a),
            "node_visits": count_nodes(left),
            "unique_nodes": len(left_unique),
        },
        "design_b": {
            "skin_thickness": float(args.skin_b),
            "spar_cap_thickness": float(args.cap_b),
            "shear_web_thickness": float(args.web_b),
            "node_visits": count_nodes(right),
            "unique_nodes": len(right_unique),
        },
        "whole_section": census,
    }
    json_path = out_dir / "csg_tree_diff.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
