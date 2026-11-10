"""Fill plot of section SDFs on a sampling grid."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from implicit_geometry_kernel.section.properties.eval_cache import SectionEvalCache
from implicit_geometry_kernel.section.boundary.labels import is_core_label
from implicit_geometry_kernel.integration.quadtree.quadrature import label_boundary_segments
from .plot_section_axes import (
    _DEFAULT_COLOR,
    _DEFAULT_CORE_ALPHA,
    _component_color,
    _naca_xy_to_plot_hv,
    _resolve_section_twist_chord,
    apply_structural_section_axes,
    blade_plot_sampling_grid,
    fit_square_section_axes,
    rotated_naca_grid_to_section_plot_axes,
    section_query_from_blade,
)
from .plot_section_hybrid import (
    plot_section_from_mesh,
    plot_section_hybrid,
)


def plot_section(section_geometry, grid, ax=None, alpha=0.45,
                 show_airfoil=False, title="Blade section geometry",
                 show_legend=True, eval_cache: SectionEvalCache | None = None,
                 *,
                 auto_limits: bool = False,
                 limit_pad_frac: float = 0.06,
                 linewidths: float = 1.0,
                 core_alpha: float | None = None,
                 omit_cores: bool = False,
                 axes_mode: str = "structural",
                 chord: float | None = None):
    """Plot all subcomponents as filled colour patches on one axis.

    Parameters
    ----------
    section_geometry : MultiCellSection or SectionMidlineHost
    grid : SDFGrid
    ax : matplotlib Axes or None
    alpha : float
        Fill transparency for structural components.
    show_airfoil : bool
        If True, overlay the chord-frame ``section_geometry.airfoil`` zero contour
        in black. Default False: that line is often misleading when the section
        is twist-rotated in body axes while ``airfoil`` stays chord-aligned.
    title : str
    show_legend : bool
        If True (default), draw a component colour legend on this axis.
    auto_limits : bool
        If True, zoom axes to the airfoil vertex bbox (plus ``limit_pad_frac``).
    limit_pad_frac : float
        Fractional pad of max(airfoil width, height) used when ``auto_limits``.
    linewidths : float
        Boundary contour width.
    core_alpha : float or None
        Override fill alpha for ``sandwich_core_*`` / ``core_*`` labels
        (defaults to ``min(0.22, alpha)``).
    omit_cores : bool
        If True, skip ``core_*`` fills so skin/spar/web read more clearly.
    axes_mode : {"structural", "default", "YB_XB_over_c"}
        ``structural`` (default) draws in canonical ``(Y_B, X_B)`` (pitch pivot
        at the origin) with the horizontal axis inverted so the LE is on the
        RHS — see ``coordinate_frames`` README. At zero twist the provided
        NACA ``grid`` is pivot-shifted; at nonzero twist a rectilinear B-frame
        sample mesh is built and SDFs are queried in rotated-NACA (remapping a
        chord-aligned lattice through twist warps cells and clips fills).
        ``default`` is an alias of ``structural`` (legacy name).
        ``YB_XB_over_c`` is the same orientation, chord-normalised.
    chord : float or None
        Chord [m] for the NACA→B pivot shift. Defaults to
        ``section_geometry.airfoil.chord``; required when that is unavailable
        or when ``axes_mode="YB_XB_over_c"``.

    Returns
    -------
    fig, ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5))
    else:
        fig = ax.figure

    mode = str(axes_mode or "structural").strip().lower()
    use_yb_xb = mode in {"yb_xb_over_c", "yb/xb/c", "scale_over_c"}
    airfoil = getattr(section_geometry, "airfoil", None)
    c = float(chord) if chord is not None else float(
        getattr(airfoil, "chord", 0.0) or 0.0
    )
    if not np.isfinite(c) or c <= 0.0:
        raise ValueError(
            "plot_section requires a positive chord (pass chord=... or set airfoil.chord)"
        )
    twist = float(getattr(section_geometry, "twist_rad", 0.0) or 0.0)
    # Zero twist: evaluate on the provided NACA grid and pivot-shift to B.
    # Nonzero twist: sample a rectilinear B mesh and query rotated-NACA SDFs —
    # remapping a chord-aligned grid through LE-rotation warps the lattice and
    # clips contourf / pcolormesh fills.
    resample_b = abs(twist) >= 1e-12
    QX = QY = None
    if resample_b:
        if airfoil is None:
            raise ValueError(
                "plot_section with nonzero twist requires section_geometry.airfoil"
            )
        ny_g, nx_g = int(np.asarray(grid.X).shape[0]), int(np.asarray(grid.X).shape[1])
        # Returns (Y_B, X_B, QX, QY); bind to X_plot/Y_plot = (horizontal, vertical).
        X_plot, Y_plot, QX, QY = blade_plot_sampling_grid(
            airfoil,
            twist_rad=twist,
            chord_m=c,
            nx=max(nx_g, 2),
            ny=max(ny_g, 2),
        )
    else:
        X_plot, Y_plot = rotated_naca_grid_to_section_plot_axes(
            grid.X, grid.Y, twist_rad=0.0, chord_m=c
        )
    if use_yb_xb:
        X_plot = np.asarray(X_plot, dtype=float) / c
        Y_plot = np.asarray(Y_plot, dtype=float) / c
        xlabel = r"$Y_B/c$ [-]"
        ylabel = r"$X_B/c$ [-]"
    else:
        xlabel = r"$Y_B$ (edgewise) [m]"
        ylabel = r"$X_B$ (flapwise) [m]"

    def _label_str(label) -> str:
        return str(getattr(label, "value", label))

    def _is_core(label) -> bool:
        # Real bay labels are ``sandwich_core_*`` (see constructive.is_core_label).
        return is_core_label(label) or _label_str(label).startswith("core")

    def _draw_priority(label) -> int:
        s = _label_str(label)
        if _is_core(label):
            return 0
        if "skin" in s:
            return 1
        if "spar" in s:
            return 2
        if "shear_web" in s or s.startswith("web"):
            return 3
        return 4

    labels = sorted(list(section_geometry), key=_draw_priority)
    if omit_cores:
        labels = [lbl for lbl in labels if not _is_core(lbl)]

    core_a = float(core_alpha) if core_alpha is not None else min(_DEFAULT_CORE_ALPHA, float(alpha))
    drawn = []
    for label in labels:
        if resample_b:
            phi = np.asarray(section_geometry[label](QX, QY), dtype=float)
        elif eval_cache is None:
            phi = grid.eval(section_geometry[label])
        else:
            phi = eval_cache.get_or_eval_with_owner(
                label, section_geometry[label], grid, owner=section_geometry
            )
        color = _component_color(label)
        fill_alpha = core_a if _is_core(label) else float(alpha)
        lw = float(linewidths) * (
            1.6 if ("shear_web" in _label_str(label) or _label_str(label).startswith("web")) else 1.0
        )
        # Filled region where phi < 0
        ax.contourf(X_plot, Y_plot, phi,
                    levels=[-np.inf, 0.0],
                    colors=[color], alpha=fill_alpha)
        # Boundary contour
        ax.contour(X_plot, Y_plot, phi,
                   levels=[0.0], colors=[color], linewidths=lw)
        drawn.append(label)

    if show_airfoil:
        if resample_b:
            phi_af = np.asarray(section_geometry.airfoil(QX, QY), dtype=float)
        elif eval_cache is None:
            phi_af = grid.eval(section_geometry.airfoil)
        else:
            phi_af = eval_cache.get_or_eval_with_owner(
                "_airfoil", section_geometry.airfoil, grid, owner=None
            )
        ax.contour(X_plot, Y_plot, phi_af,
                   levels=[0.0], colors=["k"], linewidths=2.0)

    if show_legend:
        from matplotlib.patches import Patch
        handles = [
            Patch(
                facecolor=_component_color(lbl),
                alpha=(core_a if _is_core(lbl) else float(alpha)),
                label=_label_str(lbl),
            )
            for lbl in drawn
        ]
        ax.legend(handles=handles, loc="upper right", fontsize=8)

    if auto_limits:
        # Prefer outer-skin occupancy (handles twist); fall back to airfoil verts.
        x0 = y0 = x1 = y1 = None
        skin_labels = [lbl for lbl in drawn if "skin" in _label_str(lbl)]
        limit_labels = skin_labels or list(drawn)
        for label in limit_labels:
            if resample_b:
                phi_lim = np.asarray(section_geometry[label](QX, QY), dtype=float)
            elif eval_cache is None:
                phi_lim = grid.eval(section_geometry[label])
            else:
                phi_lim = eval_cache.get_or_eval_with_owner(
                    label, section_geometry[label], grid, owner=section_geometry
                )
            mask = np.asarray(phi_lim) < 0.0
            if not np.any(mask):
                continue
            xs = np.asarray(X_plot)[mask]
            ys = np.asarray(Y_plot)[mask]
            x0, x1 = float(xs.min()), float(xs.max())
            y0, y1 = float(ys.min()), float(ys.max())
            break
        if x0 is None and getattr(section_geometry, "airfoil", None) is not None:
            verts = np.asarray(section_geometry.airfoil.vertices, dtype=float)
            if verts.ndim == 2 and verts.shape[0] >= 3:
                from implicit_geometry_kernel.frames import naca_to_twisted_pitch_axis_blade_frame
                from .frames import section_plot_axes

                h, v = section_plot_axes(
                    naca_to_twisted_pitch_axis_blade_frame(verts, twist, c)
                )
                if use_yb_xb:
                    h = h / c
                    v = v / c
                x0, x1 = float(h.min()), float(h.max())
                y0, y1 = float(v.min()), float(v.max())
        if x0 is not None:
            span = max(x1 - x0, y1 - y0, 1e-9)
            pad = float(limit_pad_frac) * span
            ax.set_xlim(x0 - pad, x1 + pad)
            ax.set_ylim(y0 - pad, y1 + pad)

    ax.set_title(title)
    apply_structural_section_axes(ax, xlabel=xlabel, ylabel=ylabel)
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

__all__ = [
    "_DEFAULT_COLOR",
    "_component_color",
    "apply_structural_section_axes",
    "blade_plot_sampling_grid",
    "fit_square_section_axes",
    "plot_section",
    "plot_section_from_mesh",
    "plot_section_hybrid",
    "plot_section_with_quadtree",
    "rotated_naca_grid_to_section_plot_axes",
    "section_query_from_blade",
]
