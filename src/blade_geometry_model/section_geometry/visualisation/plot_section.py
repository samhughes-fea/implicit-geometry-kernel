"""Section-geometry Matplotlib plotters (grid, mesh, hybrid, quadtree)."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from blade_geometry_model.section_geometry.sdf.cache import SectionEvalCache
from blade_geometry_model.area_2d.quadrature import (
    label_boundary_segments,
    label_fill_polygons,
)
from blade_geometry_model.section_geometry.constructive.labels import (
    is_core_label,
)
from blade_geometry_model.section_geometry.sections.labels import ComponentLabel


# Colour palette for subcomponents
_COMPONENT_COLORS = {
    ComponentLabel.OUTER_SKIN: "#2196F3",  # blue (legacy full annulus)
    ComponentLabel.SKIN_SUCTION_SURFACE: "#1976D2",
    ComponentLabel.SKIN_PRESSURE_SURFACE: "#64B5F6",
    ComponentLabel.SPAR_CAP_SUCTION_SURFACE: "#F44336",
    ComponentLabel.SPAR_CAP_PRESSURE_SURFACE: "#FF9800",
    "shear_web_0": "#4CAF50",
    "shear_web_1": "#009688",
    "shear_web_2": "#8BC34A",
    "sandwich_core_0": "#9C27B0",
    ComponentLabel.TRAILING_EDGE_INSERT: "#795548",
    ComponentLabel.LEADING_EDGE_INSERT: "#3F51B5",
}
_DEFAULT_COLOR = "#607D8B"
_WEB_COLORS = ["#4CAF50", "#009688", "#8BC34A", "#2E7D32"]
#: One translucent core tint for all bays (avoid dark index stacking).
_CORE_COLORS = ["#9E86C8", "#9E86C8", "#9E86C8", "#9E86C8"]
_DEFAULT_CORE_ALPHA = 0.22


def fit_square_section_axes(
    ax,
    *,
    pad_frac: float = 0.08,
    box_aspect: float | None = 1.0,
) -> None:
    """Expand current data limits to a square window (for twisted B-frame views).

    Canon owner: ``coordinate_frames.section_frames.plot.fit_square_section_axes``.
    This is a compatibility delegate.
    """
    from blade_geometry_model.coordinate_frames import (
        fit_square_section_axes as _fit_square_section_axes,
    )

    _fit_square_section_axes(ax, pad_frac=pad_frac, box_aspect=box_aspect)


def apply_structural_section_axes(
    ax,
    *,
    invert_horizontal: bool = True,
    xlabel: str | None = None,
    ylabel: str | None = None,
) -> None:
    """Finish a section plot in the canonical structural frame.

    Axes are **horizontal =** ``Y_B`` (edgewise), **vertical =** ``X_B``
    (flapwise), with the horizontal axis inverted so the LE (negative ``Y_B``)
    appears on the RHS. The inversion is display-only.

    .. important::
       **This sets up the axes; it does not transpose your data.** Which
       component to plot horizontally depends on the frame the data is in:

       * **Canonical B** (``midline_twisted_pitch_axis_blade_frame``,
         ``polyline_twisted_pitch_axis_blade_frame``,
         ``leading_edge_twisted_pitch_axis_blade_frame``) is
         stored ``(X_B, Y_B)`` — plot ``pts[:, 1]`` horizontally. Use
         ``coordinate_frames.section_plot_axes`` rather than indexing by hand.
       * **SDF grids / rotated-NACA solids** are authored with LE at
         ``x_chord = 0``. Convert to B with
         :func:`rotated_naca_grid_to_section_plot_axes` (pivot shift + twist)
         before labeling axes as ``Y_B`` / ``X_B``.

       Do **not** overlay raw NACA polylines on a B-labeled plot without
       ``coordinate_frames.naca_to_twisted_pitch_axis_blade_frame``.

    See ``geometry_utilities/coordinate_frames/README.md`` (source of truth).

    Idempotent w.r.t. horizontal inversion: safe to call after ``set_xlim``.
    """
    from blade_geometry_model.coordinate_frames import (
        setup_section_axes,
    )

    setup_section_axes(
        ax,
        xlabel=xlabel,
        ylabel=ylabel,
        invert_for_le_on_rhs=bool(invert_horizontal),
    )
    if not invert_horizontal and ax.xaxis_inverted():
        ax.invert_xaxis()


def rotated_naca_grid_to_section_plot_axes(
    X_naca: np.ndarray,
    Y_naca: np.ndarray,
    *,
    twist_rad: float = 0.0,
    chord_m: float,
    pivot_fraction: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Map SDF-grid ``(x_chord, y_thickness)`` → plot ``(Y_B, X_B)``.

    Preferred for **zero twist** (pure pitch-pivot shift). For nonzero twist,
    prefer :func:`blade_plot_sampling_grid` + SDF query in B — remapping a
    chord-aligned NACA mesh into B warps the sample lattice and clips fills.

    At zero twist: ``Y_B = x_chord - f_pivot * c``, ``X_B = y_thickness``.
    """
    from blade_geometry_model.coordinate_frames import (
        PIVOT_FRACTION_DEFAULT,
        naca_to_twisted_pitch_axis_blade_frame,
    )

    X = np.asarray(X_naca, dtype=float)
    Y = np.asarray(Y_naca, dtype=float)
    if X.shape != Y.shape:
        raise ValueError(f"X_naca and Y_naca shapes must match; got {X.shape} vs {Y.shape}")
    c_m = float(chord_m)
    if not np.isfinite(c_m) or c_m <= 0.0:
        raise ValueError(f"chord_m must be positive and finite; got {chord_m!r}")
    f_pivot = (
        PIVOT_FRACTION_DEFAULT if pivot_fraction is None else float(pivot_fraction)
    )
    t = float(twist_rad)
    qx = X.ravel()
    qy = Y.ravel()
    if abs(t) < 1e-15:
        naca = np.column_stack([qx, qy])
    else:
        # Inverse of MultiCellSection's rotate_field about the LE: R(-twist).
        ct, st = float(np.cos(t)), float(np.sin(t))
        naca = np.column_stack([ct * qx + st * qy, -st * qx + ct * qy])
    pts_b = naca_to_twisted_pitch_axis_blade_frame(naca, t, c_m, f_pivot)
    y_b = pts_b[:, 1].reshape(X.shape)
    x_b = pts_b[:, 0].reshape(X.shape)
    return y_b, x_b


def section_query_from_blade(
    pts_b: np.ndarray,
    *,
    twist_rad: float,
    chord_m: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Map canonical B ``(X_B, Y_B)`` → MultiCellSection SDF query coords.

    Twisted solids live in rotated-NACA (LE-pivot) space, not full
    ``naca_to_twisted_pitch_axis_blade_frame`` B. Convert B → NACA via
    ``twisted_pitch_axis_blade_frame_to_naca``, then apply ``R(+twist)`` about
    the LE when twisted.
    """
    from blade_geometry_model.coordinate_frames import (
        twisted_pitch_axis_blade_frame_to_naca,
    )

    naca = twisted_pitch_axis_blade_frame_to_naca(
        pts_b, float(twist_rad), float(chord_m)
    )
    t = float(twist_rad)
    if abs(t) < 1e-12:
        return naca[:, 0], naca[:, 1]
    ct, st = float(np.cos(t)), float(np.sin(t))
    qx = ct * naca[:, 0] - st * naca[:, 1]
    qy = st * naca[:, 0] + ct * naca[:, 1]
    return qx, qy


def blade_plot_sampling_grid(
    airfoil,
    *,
    twist_rad: float,
    chord_m: float | None = None,
    nx: int = 400,
    ny: int = 200,
    pad_frac: float = 0.08,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build a rectilinear B-frame sample mesh for SDF contour / pcolor plots.

    Returns
    -------
    Y_plot, X_plot, QX, QY
        ``Y_plot``/``X_plot`` are meshgrids for ``contourf(Y_plot, X_plot, …)``
        (horizontal = ``Y_B``, vertical = ``X_B``). ``QX``/``QY`` are the matching
        rotated-NACA query coordinates for ``section[label](QX, QY)``.
    """
    from blade_geometry_model.coordinate_frames import (
        naca_to_twisted_pitch_axis_blade_frame,
    )

    c_m = float(chord_m) if chord_m is not None else float(
        getattr(airfoil, "chord", 0.0) or 0.0
    )
    if not np.isfinite(c_m) or c_m <= 0.0:
        raise ValueError("blade_plot_sampling_grid requires a positive chord")
    verts = np.asarray(airfoil.vertices, dtype=float)
    verts_b = naca_to_twisted_pitch_axis_blade_frame(verts, float(twist_rad), c_m)
    span = float(max(np.ptp(verts_b[:, 0]), np.ptp(verts_b[:, 1]), 1e-9))
    pad = float(pad_frac) * span
    xb = np.linspace(float(verts_b[:, 0].min()) - pad, float(verts_b[:, 0].max()) + pad, int(ny))
    yb = np.linspace(float(verts_b[:, 1].min()) - pad, float(verts_b[:, 1].max()) + pad, int(nx))
    # indexing='xy': Y_plot[i,j] varies with j (horizontal), X_plot with i (vertical)
    Y_plot, X_plot = np.meshgrid(yb, xb)
    pts_b = np.column_stack([X_plot.ravel(), Y_plot.ravel()])
    qx, qy = section_query_from_blade(
        pts_b, twist_rad=float(twist_rad), chord_m=c_m
    )
    return Y_plot, X_plot, qx.reshape(X_plot.shape), qy.reshape(X_plot.shape)


def _resolve_section_twist_chord(
    section=None,
    mesh=None,
    *,
    twist_rad: float | None = None,
    chord_m: float | None = None,
) -> tuple[float, float]:
    """Resolve twist [rad] and chord [m] for NACA→B plot remaps."""
    twist = float(twist_rad) if twist_rad is not None else float(
        getattr(section, "twist_rad", 0.0) or 0.0
    )
    c = float(chord_m) if chord_m is not None else 0.0
    if not (np.isfinite(c) and c > 0.0) and section is not None:
        af = getattr(section, "airfoil", None) or getattr(section, "_af", None)
        c = float(getattr(af, "chord", 0.0) or 0.0)
    if not (np.isfinite(c) and c > 0.0) and mesh is not None:
        spec = getattr(mesh, "spec", None)
        if spec is not None:
            # Quadtree domain is NACA-like; chord ≈ LE→TE span on x.
            c = float(getattr(spec, "x_max", 0.0) - getattr(spec, "x_min", 0.0))
    if not (np.isfinite(c) and c > 0.0):
        raise ValueError(
            "chord_m is required for B-frame section plots "
            "(pass chord_m=... or attach airfoil.chord)"
        )
    return twist, c


def _naca_xy_to_plot_hv(
    xy: np.ndarray,
    *,
    twist_rad: float,
    chord_m: float,
) -> np.ndarray:
    """``(N, 2)`` rotated-NACA → plot coords ``(N, 2)`` as ``(Y_B, X_B)``."""
    p = np.asarray(xy, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2:
        raise ValueError(f"xy must be (N, 2); got {p.shape}")
    h, v = rotated_naca_grid_to_section_plot_axes(
        p[:, 0], p[:, 1], twist_rad=float(twist_rad), chord_m=float(chord_m)
    )
    return np.column_stack([np.asarray(h, dtype=float), np.asarray(v, dtype=float)])


def _component_color(label):
    key = str(getattr(label, "value", label))
    if label in _COMPONENT_COLORS:
        return _COMPONENT_COLORS[label]
    if key in _COMPONENT_COLORS:
        return _COMPONENT_COLORS[key]
    if key.startswith("skin_suction_surface") or key == "skin_suction_surface":
        return _COMPONENT_COLORS[ComponentLabel.SKIN_SUCTION_SURFACE]
    if key.startswith("skin_pressure_surface") or key == "skin_pressure_surface":
        return _COMPONENT_COLORS[ComponentLabel.SKIN_PRESSURE_SURFACE]
    if key == "outer_skin" or key.startswith("outer_skin_"):
        return _COMPONENT_COLORS[ComponentLabel.OUTER_SKIN]
    if key.startswith("shear_web_"):
        try:
            idx = int(key.rsplit("_", 1)[-1])
            return _WEB_COLORS[idx % len(_WEB_COLORS)]
        except (ValueError, IndexError):
            return _WEB_COLORS[0]
    if key.startswith("sandwich_core_"):
        suffix = key.rsplit("_", 1)[-1]
        try:
            idx = int(suffix)
            return _CORE_COLORS[idx % len(_CORE_COLORS)]
        except ValueError:
            return _CORE_COLORS[0]
    return _COMPONENT_COLORS.get(label, _COMPONENT_COLORS.get(key, _DEFAULT_COLOR))


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
                from blade_geometry_model.coordinate_frames import (
                    naca_to_twisted_pitch_axis_blade_frame,
                    section_plot_axes,
                )

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
    """Plot wall solids from constructive polygons and residuals from quadtree mesh.

    Constructive rings and residual mesh leaves are rotated-NACA; with
    ``structural_axes`` they are remapped to canonical B plot axes.

    Sandwich cores use ``core_alpha`` (default ``min(0.22, alpha)``) and are
    drawn under wall solids so translucent bay fills do not cover skins/webs.
    """
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.patches import Patch, PathPatch
    from matplotlib.path import Path

    from blade_geometry_model.section_geometry.constructive import (
        emit_constructive_polygons,
        is_residual_label,
    )

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
        from shapely.geometry.polygon import orient

        g = orient(g, sign=1.0)
        verts_list = []
        codes_list = []
        for ring in (g.exterior, *g.interiors):
            c = _to_plot(np.asarray(ring.coords, dtype=float))
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
