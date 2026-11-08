"""Mesh, hybrid polygon, and quadtree-panel section plotters."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from implicit_geometry_kernel.integration.quadtree.quadrature import (
    label_boundary_segments,
    label_fill_polygons,
)
from implicit_geometry_kernel.section.boundary.labels import is_core_label
from .plot_section_axes import (
    _DEFAULT_CORE_ALPHA,
    _component_color,
    _naca_xy_to_plot_hv,
    _resolve_section_twist_chord,
    apply_structural_section_axes,
)


def plot_section_from_mesh(mesh, labels=None, ax=None, *, alpha=0.45,
                           title="Blade section geometry", show_legend=True,
                           color_fn=_component_color, airfoil_vertices=None,
                           structural_axes: bool = True,
                           twist_rad: float | None = None,
                           chord_m: float | None = None,
                           section=None):
    """Plot section components from a shared quadtree mesh.

    Mesh coordinates are rotated-NACA (LE origin). With ``structural_axes`` they
    are remapped to canonical B plot axes (pitch pivot at the origin).
    """
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.patches import Patch

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5))
    else:
        fig = ax.figure

    if labels is None:
        labels = list(mesh.phi_by_label.keys())
    labels = [str(getattr(lbl, "value", lbl)) for lbl in labels]

    twist = 0.0
    chord = 1.0
    if structural_axes:
        twist, chord = _resolve_section_twist_chord(
            section, mesh, twist_rad=twist_rad, chord_m=chord_m
        )

    def _to_plot(poly: np.ndarray) -> np.ndarray:
        if not structural_axes:
            return np.asarray(poly, dtype=float)
        return _naca_xy_to_plot_hv(poly, twist_rad=twist, chord_m=chord)

    for label in labels:
        color = color_fn(label)
        polys = [_to_plot(p) for p in label_fill_polygons(mesh, label)]
        if polys:
            ax.add_collection(
                PolyCollection(
                    polys, facecolors=color, edgecolors="none",
                    alpha=alpha, rasterized=True,
                )
            )
        segs = label_boundary_segments(mesh, label)
        if len(segs):
            segs_p = np.asarray([_to_plot(s) for s in segs], dtype=float)
            ax.add_collection(LineCollection(segs_p, colors=color, linewidths=1.0))

    if airfoil_vertices is not None:
        v = np.asarray(airfoil_vertices, dtype=float)
        if v.ndim == 2 and v.shape[0] >= 2:
            vv = _to_plot(np.vstack([v, v[:1]]))
            ax.plot(vv[:, 0], vv[:, 1], color="k", linewidth=2.0)

    if show_legend:
        handles = [Patch(facecolor=color_fn(lbl), alpha=alpha, label=lbl) for lbl in labels]
        ax.legend(handles=handles, loc="upper right", fontsize=8)

    spec = mesh.spec
    corners = np.array(
        [
            [spec.x_min, spec.y_min],
            [spec.x_max, spec.y_min],
            [spec.x_max, spec.y_max],
            [spec.x_min, spec.y_max],
        ],
        dtype=float,
    )
    corners_p = _to_plot(corners)
    ax.set_xlim(float(corners_p[:, 0].min()), float(corners_p[:, 0].max()))
    ax.set_ylim(float(corners_p[:, 1].min()), float(corners_p[:, 1].max()))
    ax.set_title(title)
    if structural_axes:
        apply_structural_section_axes(ax)
    else:
        ax.set_aspect("equal")
        ax.set_xlabel("x_chord [m]")
        ax.set_ylabel("y_thickness [m]")
    return fig, ax


def plot_section_hybrid(
    section,
    residual_mesh=None,
    ax=None,
    *,
    alpha=0.45,
    core_alpha: float | None = None,
    omit_cores: bool = False,
    title="Blade section geometry (hybrid)",
    show_legend=True,
    color_fn=_component_color,
    structural_axes: bool = True,
    twist_rad: float | None = None,
    chord_m: float | None = None,
):
    """Wall solids from constructive polygons; residuals from the quadtree mesh."""
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.patches import Patch, PathPatch
    from matplotlib.path import Path

    from implicit_geometry_kernel.section.boundary.cores import emit_constructive_polygons

    from implicit_geometry_kernel.section.boundary.labels import is_residual_label

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5))
    else:
        fig = ax.figure

    twist = 0.0
    chord = 1.0
    if structural_axes:
        twist, chord = _resolve_section_twist_chord(
            section, residual_mesh, twist_rad=twist_rad, chord_m=chord_m
        )

    def _to_plot(xy: np.ndarray) -> np.ndarray:
        if not structural_axes:
            return np.asarray(xy, dtype=float)
        return _naca_xy_to_plot_hv(xy, twist_rad=twist, chord_m=chord)

    def _polygon_path(g) -> Path:
        from implicit_geometry_kernel.section.boundary.polygon_ops import oriented_rings

        verts_list = []
        codes_list = []
        for ring in oriented_rings(g):
            c = _to_plot(np.asarray(ring, dtype=float))
            verts_list.append(c)
            codes = np.full(len(c), Path.LINETO, dtype=Path.code_type)
            codes[0] = Path.MOVETO
            codes[-1] = Path.CLOSEPOLY
            codes_list.append(codes)
        return Path(np.vstack(verts_list), np.concatenate(codes_list))

    def _draw_priority(label: str) -> int:
        if is_core_label(label):
            return 0
        if "skin" in label:
            return 1
        if "spar" in label:
            return 2
        if "shear_web" in label or label.startswith("web"):
            return 3
        return 4

    core_a = float(core_alpha) if core_alpha is not None else min(_DEFAULT_CORE_ALPHA, float(alpha))
    legend_labels: list[str] = []
    polys = emit_constructive_polygons(section)
    if omit_cores:
        polys = {k: v for k, v in polys.items() if not is_core_label(k)}
    ordered = sorted(polys.items(), key=lambda kv: _draw_priority(str(kv[0])))
    xs_lim: list[float] = []
    ys_lim: list[float] = []
    for label, geom in ordered:
        lab = str(label)
        color = color_fn(label)
        legend_labels.append(lab)
        fill_a = core_a if is_core_label(lab) else float(alpha)
        z = 1 if is_core_label(lab) else 2
        geoms = getattr(geom, "geoms", None) or (geom,)
        for g in geoms:
            if g.is_empty:
                continue
            path = _polygon_path(g)
            xs_lim.extend(path.vertices[:, 0].tolist())
            ys_lim.extend(path.vertices[:, 1].tolist())
            ax.add_patch(
                PathPatch(
                    path,
                    facecolor=color,
                    edgecolor="none",
                    alpha=fill_a,
                    rasterized=True,
                    zorder=z,
                )
            )

    if residual_mesh is not None and getattr(residual_mesh, "phi_by_label", None):
        for label in residual_mesh.phi_by_label:
            lab = str(getattr(label, "value", label))
            if not is_residual_label(lab):
                continue
            color = color_fn(lab)
            legend_labels.append(lab)
            fills = [_to_plot(p) for p in label_fill_polygons(residual_mesh, label)]
            if fills:
                ax.add_collection(
                    PolyCollection(
                        fills,
                        facecolors=color,
                        edgecolors="none",
                        alpha=alpha,
                        rasterized=True,
                        zorder=2,
                    )
                )
            segs = label_boundary_segments(residual_mesh, label)
            if len(segs):
                segs_p = np.asarray([_to_plot(s) for s in segs], dtype=float)
                ax.add_collection(LineCollection(segs_p, colors=color, linewidths=1.0, zorder=3))
        spec = residual_mesh.spec
        corners = _to_plot(
            np.array(
                [
                    [spec.x_min, spec.y_min],
                    [spec.x_max, spec.y_min],
                    [spec.x_max, spec.y_max],
                    [spec.x_min, spec.y_max],
                ],
                dtype=float,
            )
        )
        xs_lim.extend(corners[:, 0].tolist())
        ys_lim.extend(corners[:, 1].tolist())

    if xs_lim and ys_lim:
        pad = 0.05 * max(
            max(xs_lim) - min(xs_lim),
            max(ys_lim) - min(ys_lim),
            1e-6,
        )
        ax.set_xlim(min(xs_lim) - pad, max(xs_lim) + pad)
        ax.set_ylim(min(ys_lim) - pad, max(ys_lim) + pad)
    else:
        host = getattr(section, "airfoil", None) or getattr(section, "_af", None)
        verts = getattr(host, "vertices", None)
        if verts is not None:
            v = _to_plot(np.asarray(verts, dtype=float))
            pad = 0.05 * max(np.ptp(v[:, 0]), np.ptp(v[:, 1]), 1e-6)
            ax.set_xlim(float(v[:, 0].min()) - pad, float(v[:, 0].max()) + pad)
            ax.set_ylim(float(v[:, 1].min()) - pad, float(v[:, 1].max()) + pad)

    if show_legend and legend_labels:
        handles = [
            Patch(
                facecolor=color_fn(lbl),
                alpha=(core_a if is_core_label(lbl) else float(alpha)),
                label=lbl,
            )
            for lbl in legend_labels
        ]
        ax.legend(handles=handles, loc="upper right", fontsize=8)

    ax.set_title(title)
    if structural_axes:
        apply_structural_section_axes(ax)
    else:
        ax.set_aspect("equal")
        ax.set_xlabel("x_chord [m]")
        ax.set_ylabel("y_thickness [m]")
    return fig, ax


def plot_section_with_quadtree(
    mesh,
    labels=None,
    *,
    title="Blade section geometry",
    mesh_title=None,
    alpha=0.45,
    show_legend=True,
    color_fn=_component_color,
    airfoil_vertices=None,
    figsize=(20, 6),
    section=None,
    twist_rad: float | None = None,
    chord_m: float | None = None,
):
    """Two-panel figure: built geometry (left) and union-refined quadtree mesh (right)."""
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Polygon

    if labels is None:
        labels = list(mesh.phi_by_label.keys())
    labels = [str(getattr(lbl, "value", lbl)) for lbl in labels]

    twist, chord = _resolve_section_twist_chord(
        section, mesh, twist_rad=twist_rad, chord_m=chord_m
    )

    leaves = mesh.leaf_cells()
    if mesh_title is None:
        max_depth = max((int(getattr(c, "level", 0)) for c in leaves), default=0)
        mesh_title = (
            f"Square union-refined quadtree (depth<={max_depth}, leaves={len(leaves)})"
        )

    fig, (ax_geom, ax_quad) = plt.subplots(1, 2, figsize=figsize)
    if section is not None:
        plot_section_hybrid(
            section,
            residual_mesh=mesh,
            ax=ax_geom,
            alpha=alpha,
            title=title,
            show_legend=show_legend,
            color_fn=color_fn,
            twist_rad=twist,
            chord_m=chord,
        )
    else:
        plot_section_from_mesh(
            mesh,
            labels=labels,
            ax=ax_geom,
            alpha=alpha,
            title=title,
            show_legend=show_legend,
            color_fn=color_fn,
            airfoil_vertices=airfoil_vertices,
            twist_rad=twist,
            chord_m=chord,
            section=section,
        )

    for cell in leaves:
        corners = _naca_xy_to_plot_hv(
            np.array(
                [
                    [cell.x_min, cell.y_min],
                    [cell.x_max, cell.y_min],
                    [cell.x_max, cell.y_max],
                    [cell.x_min, cell.y_max],
                ],
                dtype=float,
            ),
            twist_rad=twist,
            chord_m=chord,
        )
        ax_quad.add_patch(
            Polygon(
                corners,
                fill=False,
                edgecolor="0.6",
                linewidth=0.3,
                closed=True,
            )
        )
    for label in labels:
        segs = label_boundary_segments(mesh, label)
        if len(segs):
            segs_p = np.asarray(
                [_naca_xy_to_plot_hv(s, twist_rad=twist, chord_m=chord) for s in segs],
                dtype=float,
            )
            ax_quad.add_collection(
                LineCollection(segs_p, colors=[color_fn(label)], linewidths=1.4)
            )

    spec = mesh.spec
    corners = _naca_xy_to_plot_hv(
        np.array(
            [
                [spec.x_min, spec.y_min],
                [spec.x_max, spec.y_min],
                [spec.x_max, spec.y_max],
                [spec.x_min, spec.y_max],
            ],
            dtype=float,
        ),
        twist_rad=twist,
        chord_m=chord,
    )
    ax_quad.set_xlim(float(corners[:, 0].min()), float(corners[:, 0].max()))
    ax_quad.set_ylim(float(corners[:, 1].min()), float(corners[:, 1].max()))
    ax_quad.set_title(mesh_title)
    apply_structural_section_axes(ax_quad)
    if section is not None:
        ax_geom.set_xlim(ax_quad.get_xlim())
        ax_geom.set_ylim(ax_quad.get_ylim())
        apply_structural_section_axes(ax_geom)
    fig.tight_layout()
    return fig, (ax_geom, ax_quad)


__all__ = ["plot_section_from_mesh", "plot_section_hybrid", "plot_section_with_quadtree"]
