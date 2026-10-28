"""
sections.multicell
==================
MultiCellSection — generalised N-web blade section with N+1 sandwich-core bays.
"""

from __future__ import annotations

from blade_geometry_model.section_contract.errors import (
    InfeasibleSectionGeometry,
)

from typing import Any

import numpy as np

from blade_geometry_model.section_geometry.sdf.ir import Offset, Rotate, union_expr
from blade_geometry_model.section_geometry.sdf.profiling import timed_stage
from blade_geometry_model.section_geometry.sdf.transforms import rotate_field

from ..layout.thickness_limits import (
    clamp_skin_thickness_m,
    clamp_spar_thickness_m,
    clamp_web_thickness_m,
)
from ..layout.cap_placement import CapPlacementRule
from ..layout.structural import parse_fixed_cap_anchor, parse_structural_family
from .clip_algebra import (
    _parallel_web_half_gt,
    _parallel_web_half_lt,
    _parallel_web_strip_clip,
    airfoil_expr,
    flapwise_web_normal,
    web_station_projection,
)
from blade_geometry_model.section_geometry.recipes.sdf_ir.base import (
    expr_callable,
    twist_component,
)
from .build_spar_caps import build_spar_caps
from .labels import ComponentLabel
from .web_layout import flapwise_web_sframe_inputs, resolve_web_anchors


def _label_str(label: Any) -> str:
    return str(getattr(label, "value", label))


class MultiCellSection:
    """N-web generalised blade section (N+1 structural bays, N+1 foam cores when enabled)."""

    def __init__(
        self,
        airfoil_sdf,
        web_x_positions,
        web_thickness=0.004,
        web_alignment="chord_normal",
        cap_height=0.012,
        skin_thickness=0.003,
        twist_angle=0.0,
        web_y_coords=None,
        trailing_edge_insert_x_m=None,
        leading_edge_insert_x_m=None,
        le_radius=None,
        core_enabled=True,
        y_search=None,
        structural_family="D",
        fixed_cap_anchor="pitching",
        pitch_fraction_of_chord_from_le=1.0 / 3.0,
        fixed_cap_chord_half_width=None,
        discrete_cap_chord_half_width=None,
        *,
        cap_placement=None,
        web_endpoints=None,
        stage_sink: dict[str, float] | None = None,
        spar_cap_ply_pitch_m: float | None = None,
    ):
        """Assemble a multi-cell section.

        ``web_endpoints`` (optional): sequence of ``(x_top, y_top, x_bot, y_bot)``
        for non-parallel / slanted webs. When set, ``web_x_positions`` /
        ``web_y_coords`` are ignored for axis geometry (bay clips still use the
        chordwise mid-station of each axis). True Y-branching (three or more
        walls at one node) is out of scope.

        ``spar_cap_ply_pitch_m`` (optional): floor post-clamp cap heights onto
        this ply pitch so CapWall / midlines stay on-grid for Catalog FI bind
        (BL-019).
        """
        explicit_endpoints = None
        if web_endpoints is not None:
            explicit_endpoints = [
                (float(xt), float(yt), float(xb), float(yb))
                for xt, yt, xb, yb in web_endpoints
            ]
            xs_requested = [
                0.5 * (e[0] + e[2]) for e in explicit_endpoints
            ]
        else:
            xs_requested = sorted(float(x) for x in web_x_positions)
        # Catalogue echo only — assemble control flow uses CapPlacementRule.
        self._structural_family = parse_structural_family(structural_family)
        if cap_placement is not None:
            _cap_rule = cap_placement
        else:
            _pos = None
            fam = self._structural_family
            if fam == "B":
                _pos = parse_fixed_cap_anchor(fixed_cap_anchor)
            _cap_rule = CapPlacementRule.from_structural_family(
                fam, spar_cap_positioning=_pos
            )
        if len(xs_requested) < 1 and _cap_rule.kind not in ("none", "fixed_band"):
            raise ValueError(
                "web_x_positions is empty, which is only supported for "
                f"cap placement 'none' or 'fixed_band' (got {_cap_rule.kind!r}). "
                "per_web_pads needs one pad per web; continuous_box spans "
                "outermost webs."
            )

        skin_arg = skin_thickness
        if isinstance(skin_arg, (int, float)):
            t_ss = t_ps = clamp_skin_thickness_m(skin_arg)
        else:
            t_ss = clamp_skin_thickness_m(float(skin_arg[0]))
            t_ps = clamp_skin_thickness_m(float(skin_arg[1]))
        skin_thickness = max(t_ss, t_ps)

        self._af = airfoil_sdf
        self._skin_t = float(skin_thickness)
        self._skin_t_ss = float(t_ss)
        self._skin_t_ps = float(t_ps)
        self._twist = float(twist_angle)
        if _cap_rule.kind in ("fixed_band", "per_web_pads"):
            from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF as _AirfoilSDF

            if not isinstance(airfoil_sdf, _AirfoilSDF):
                raise TypeError(
                    "fixed_band / per_web_pads require airfoil_sdf "
                    "to be AirfoilSDF (pitching / max-thickness anchors)."
                )

        self._airfoil_expr = airfoil_expr(airfoil_sdf)
        if self._airfoil_expr is None:
            raise ValueError("MultiCellSection requires AirfoilSDF with vertices.")

        self._components: dict[str, Any] = {}
        self._cap_placement = _cap_rule

        _chord = float(getattr(airfoil_sdf, "chord", 1.0))
        _pitch_fr = float(pitch_fraction_of_chord_from_le)
        _fhw = (
            float(fixed_cap_chord_half_width)
            if fixed_cap_chord_half_width is not None
            else (
                float(_cap_rule.half_width)
                if _cap_rule.half_width is not None
                else 0.05 * _chord
            )
        )
        if _cap_rule.kind == "fixed_band" and _cap_rule.half_width is None:
            _cap_rule = CapPlacementRule.fixed_band(
                anchor=_cap_rule.anchor or "max_thickness",
                half_width=_fhw,
            )
            self._cap_placement = _cap_rule
        _dhw = (
            float(discrete_cap_chord_half_width)
            if discrete_cap_chord_half_width is not None
            else 0.04 * _chord
        )

        with timed_stage("assemble_skin", stage_sink):
            from blade_geometry_model.section_geometry.recipes.walls import (
                SkinWall,
            )

            # SS/PS half-annuli (suction = camber-upper, pressure = camber-lower).
            wall_recipes: dict[str, Any] = {}
            skin_ss_wall = SkinWall(
                thickness=float(t_ss), skin_side="suction_surface"
            )
            skin_ps_wall = SkinWall(
                thickness=float(t_ps), skin_side="pressure_surface"
            )
            skin_ss = skin_ss_wall.to_sdf_component(
                airfoil_sdf, airfoil_expr_node=self._airfoil_expr
            )
            skin_ps = skin_ps_wall.to_sdf_component(
                airfoil_sdf, airfoil_expr_node=self._airfoil_expr
            )
            wall_recipes[_label_str(ComponentLabel.SKIN_SUCTION_SURFACE)] = skin_ss_wall
            wall_recipes[_label_str(ComponentLabel.SKIN_PRESSURE_SURFACE)] = skin_ps_wall
            self._components[ComponentLabel.SKIN_SUCTION_SURFACE] = skin_ss
            self._components[ComponentLabel.SKIN_PRESSURE_SURFACE] = skin_ps
            # Boundary SDFs stay full-annulus (cores/webs/inserts use the cavity).
            skin_outer_sdf = skin_ss.outer_sdf
            skin_inner_sdf = skin_ss.inner_sdf
            self._skin_outer_boundary_expr = Offset(self._airfoil_expr, amount=0.0)
            self._skin_inner_boundary_expr = Offset(self._airfoil_expr, amount=-skin_thickness)
            phi_inner_ref = expr_callable(Offset(self._airfoil_expr, amount=-skin_thickness))

        with timed_stage("assemble_anchors", stage_sink):
            if explicit_endpoints is not None:
                xs = [0.5 * (e[0] + e[2]) for e in explicit_endpoints]
                anchors = [(e[1], e[3]) for e in explicit_endpoints]
            elif web_y_coords is None:
                # Assemble is fail-closed: every *requested* web must anchor.
                # Partial drops belong to resolve_station_layout / freeze, which
                # step the catalogue (2D-F -> 1C-F) *before* assemble. Shrinking
                # xs here while keeping continuous_box built illegal 1D solids
                # (None caps -> twist AttributeError; Station 38).
                resolved, dropped = resolve_web_anchors(
                    airfoil_sdf,
                    skin_thickness,
                    xs_requested,
                    y_search=y_search,
                    phi_inner=phi_inner_ref,
                    return_dropped=True,
                )
                if xs_requested and not resolved:
                    raise InfeasibleSectionGeometry(
                        "No feasible web anchors: local thickness is smaller than "
                        f"inner-skin offset (skin_thickness={skin_thickness:.6g}) "
                        f"at requested x={xs_requested}."
                    )
                if dropped:
                    drop_desc = ", ".join(
                        f"{float(d.x_m):.6g} ({d.reason})" for d in dropped
                    )
                    raise InfeasibleSectionGeometry(
                        "Requested webs cannot all anchor at this skin_thickness="
                        f"{skin_thickness:.6g}: dropped [{drop_desc}]. "
                        "Call resolve_station_layout (or pass frozen web_y_coords) "
                        "so the layout steps down before assemble; assemble does "
                        "not silently shrink the web set."
                    )
                xs = [a[0] for a in resolved]
                anchors = [(a[1], a[2]) for a in resolved]
            else:
                # Frozen anchors from resolve_station_layout / sizing freeze:
                # do not re-probe in a way that can shorten the web set.
                xs = list(xs_requested)
                anchors = [(float(yt), float(yb)) for yt, yb in web_y_coords]
                if len(anchors) != len(xs):
                    raise ValueError(f"web_y_coords must have {len(xs)} entries.")

        N = len(xs)
        # Catalogue has no 1D: a continuous box needs fore and aft webs.
        if _cap_rule.kind == "continuous_box" and N < 2:
            raise InfeasibleSectionGeometry(
                f"continuous_box cap placement requires at least 2 webs (got N={N}). "
                "There is no 1D catalogue layout; resolve_station_layout must step "
                "D-family to 1C (per_web_pads) before assemble."
            )
        if isinstance(web_thickness, (int, float)):
            web_thicknesses = [float(web_thickness)] * N
        else:
            web_thicknesses = [float(t) for t in web_thickness]
            if len(web_thicknesses) != N:
                raise ValueError(
                    f"web_thickness must be scalar or length {N} (got {len(web_thicknesses)})."
                )
        web_thicknesses = [clamp_web_thickness_m(t) for t in web_thicknesses]

        if isinstance(web_alignment, str):
            web_alignments = [web_alignment] * N
        else:
            web_alignments = list(web_alignment)
            if len(web_alignments) != N:
                raise ValueError(
                    f"web_alignment must be scalar or length {N} (got {len(web_alignments)})."
                )

        if isinstance(cap_height, (int, float)):
            cap_h_upper = cap_h_lower = float(cap_height)
        else:
            cap_h_upper, cap_h_lower = float(cap_height[0]), float(cap_height[1])
        cap_h_upper = clamp_spar_thickness_m(cap_h_upper)
        cap_h_lower = clamp_spar_thickness_m(cap_h_lower)

        with timed_stage("assemble_webs", stage_sink):
            all_flapwise = N > 0 and all(a == "flapwise" for a in web_alignments)
            if all_flapwise:
                _nx, _ny = flapwise_web_normal(float(twist_angle))
                c_station = [
                    web_station_projection(
                        _nx, _ny, xs[i], 0.5 * (anchors[i][0] + anchors[i][1])
                    )
                    for i in range(N)
                ]
            else:
                _nx = _ny = None
                c_station = None

            flapwise_present = any(a == "flapwise" for a in web_alignments)
            needs_b_frame_lookup = flapwise_present and abs(float(twist_angle)) > 1e-10
            phi_inner_b = None
            y_search_b = y_search
            if needs_b_frame_lookup:
                phi_inner_b = rotate_field(phi_inner_ref, float(twist_angle))
                verts = getattr(airfoil_sdf, "vertices", None)
                if verts is not None:
                    arr = np.asarray(verts, dtype=float)
                    cos_t = float(np.cos(float(twist_angle)))
                    sin_t = float(np.sin(float(twist_angle)))
                    y_b_vals = sin_t * arr[:, 0] + cos_t * arr[:, 1]
                    y_search_b = (float(y_b_vals.min()) - 0.1, float(y_b_vals.max()) + 0.1)
                elif y_search_b is None:
                    y_search_b = (-1.5, 1.5)

            web_components = []
            from blade_geometry_model.section_geometry.recipes.walls import (
                WebWall,
            )

            for i, (x, t, align, (y_top, y_bot)) in enumerate(
                zip(xs, web_thicknesses, web_alignments, anchors)
            ):
                if explicit_endpoints is not None:
                    x_top_S, y_top_S, x_bot_S, y_bot_S = explicit_endpoints[i]
                elif align == "flapwise" and needs_b_frame_lookup:
                    x_top_S, y_top_S, x_bot_S, y_bot_S = flapwise_web_sframe_inputs(
                        x, y_top, y_bot, twist_angle, phi_inner_b, y_search_b
                    )
                else:
                    x_top_S = x
                    x_bot_S = x
                    y_top_S = y_top
                    y_bot_S = y_bot

                web_wall = WebWall(
                    skin_thickness=float(skin_thickness),
                    x_top=float(x_top_S),
                    y_top=float(y_top_S),
                    x_bot=float(x_bot_S),
                    y_bot=float(y_bot_S),
                    thickness=float(t),
                    alignment=str(align),
                    twist_angle=float(twist_angle),
                )
                label = f"shear_web_{i}"
                web = web_wall.to_sdf_component(
                    airfoil_sdf, airfoil_expr_node=self._airfoil_expr
                )
                wall_recipes[label] = web_wall
                self._components[label] = web
                web_components.append(web)

        cap_upper = cap_lower = None
        spar_cap_components_unrotated: dict[str, Any] = {}
        spar_cap_wall_recipes: dict[str, Any] = {}

        with timed_stage("assemble_caps", stage_sink):
            if _cap_rule.kind != "none":
                (
                    cap_upper,
                    cap_lower,
                    spar_cap_components_unrotated,
                    spar_cap_wall_recipes,
                ) = build_spar_caps(
                    airfoil_sdf,
                    skin_thickness,
                    xs,
                    anchors,
                    cap_h_upper,
                    cap_h_lower,
                    web_components,
                    cap_placement=_cap_rule,
                    airfoil_expr_node=self._airfoil_expr,
                    all_flapwise=all_flapwise,
                    c_station=c_station,
                    nx=_nx,
                    ny=_ny,
                    N=N,
                    pitch_fr=_pitch_fr,
                    fhw=_fhw,
                    dhw=_dhw,
                    y_search=y_search,
                    phi_inner_ref=phi_inner_ref,
                    spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
                )
                # Cap run-out (BL-006) is decided *before* assemble via
                # layout_without_spar_caps. Never store None placeholders:
                # twist_component would crash on .expr (Station 38).
                if cap_upper is None or cap_lower is None:
                    raise InfeasibleSectionGeometry(
                        "Cap placement "
                        f"{_cap_rule.kind!r} expected spar caps but build_spar_caps "
                        f"returned suction={cap_upper is not None} "
                        f"pressure={cap_lower is not None}. Use "
                        "layout_without_spar_caps when the cavity cannot hold "
                        "minimum-gauge caps (BL-006); do not assemble None caps."
                    )
                self._components[ComponentLabel.SPAR_CAP_SUCTION_SURFACE] = cap_upper
                self._components[ComponentLabel.SPAR_CAP_PRESSURE_SURFACE] = cap_lower
                for label, wall in spar_cap_wall_recipes.items():
                    wall_recipes[label] = wall

        self._spar_cap_components_unrotated = dict(spar_cap_components_unrotated)
        self._spar_cap_wall_recipes = dict(spar_cap_wall_recipes)

        # Wall solids excluded from sandwich-core fill (not layup/material).
        wall_solids_for_core = (
            list(web_components)
            if _cap_rule.kind == "none"
            else [cap_upper, cap_lower] + web_components
        )
        wall_exclusion_exprs = [c.expr for c in wall_solids_for_core]

        with timed_stage("assemble_cores", stage_sink):
            from blade_geometry_model.section_geometry.recipes.walls import (
                CoreWall,
            )

            if core_enabled and N >= 1:
                n_bays = N + 1
                for k in range(n_bays):
                    if all_flapwise and c_station is not None:
                        if k == 0:
                            sclip = _parallel_web_half_lt(_nx, _ny, c_station[0])
                            fw_mode = "half_lt"
                            c_lo, c_hi = float(c_station[0]), None
                        elif k == n_bays - 1:
                            sclip = _parallel_web_half_gt(_nx, _ny, c_station[-1])
                            fw_mode = "half_gt"
                            c_lo, c_hi = float(c_station[-1]), None
                        else:
                            sclip = _parallel_web_strip_clip(
                                _nx, _ny, c_station[k - 1], c_station[k]
                            )
                            fw_mode = "strip"
                            c_lo = float(min(c_station[k - 1], c_station[k]))
                            c_hi = float(max(c_station[k - 1], c_station[k]))
                        core_wall = CoreWall(
                            skin_thickness=float(skin_thickness),
                            bay_mode="flapwise",
                            flapwise_nx=float(_nx),
                            flapwise_ny=float(_ny),
                            flapwise_c_lo=c_lo,
                            flapwise_c_hi=c_hi,
                            flapwise_clip_mode=fw_mode,  # type: ignore[arg-type]
                        )
                        core = core_wall.to_sdf_component(
                            airfoil_sdf,
                            airfoil_expr_node=self._airfoil_expr,
                            exclusion_exprs=wall_exclusion_exprs,
                            strip_clip=sclip,
                        )
                    else:
                        if k == 0:
                            x_s, x_e = None, xs[0]
                        elif k == n_bays - 1:
                            x_s, x_e = xs[-1], None
                        else:
                            x_s, x_e = xs[k - 1], xs[k]
                        core_wall = CoreWall(
                            skin_thickness=float(skin_thickness),
                            bay_mode="chordwise",
                            x_start=x_s,
                            x_end=x_e,
                        )
                        core = core_wall.to_sdf_component(
                            airfoil_sdf,
                            airfoil_expr_node=self._airfoil_expr,
                            exclusion_exprs=wall_exclusion_exprs,
                        )
                    label = f"sandwich_core_{k}"
                    wall_recipes[label] = core_wall
                    self._components[label] = core
            elif core_enabled and N == 0:
                # Zero webs: a single whole-envelope bay (no chordwise/flapwise
                # split point exists), excluding whatever cap regions were built.
                core_wall = CoreWall(
                    skin_thickness=float(skin_thickness),
                    bay_mode="chordwise",
                    x_start=None,
                    x_end=None,
                )
                core = core_wall.to_sdf_component(
                    airfoil_sdf,
                    airfoil_expr_node=self._airfoil_expr,
                    exclusion_exprs=wall_exclusion_exprs,
                )
                wall_recipes["sandwich_core_0"] = core_wall
                self._components["sandwich_core_0"] = core

        if trailing_edge_insert_x_m is not None:
            from blade_geometry_model.section_geometry.recipes.walls import (
                InsertWall,
            )

            te_wall = InsertWall(
                kind="trailing_edge",
                skin_thickness=float(skin_thickness),
                x_limit=float(trailing_edge_insert_x_m),
            )
            te = te_wall.to_sdf_component(
                airfoil_sdf, airfoil_expr_node=self._airfoil_expr
            )
            wall_recipes[_label_str(ComponentLabel.TRAILING_EDGE_INSERT)] = te_wall
            self._components[ComponentLabel.TRAILING_EDGE_INSERT] = te

        if leading_edge_insert_x_m is not None:
            from blade_geometry_model.section_geometry.recipes.walls import (
                InsertWall,
            )

            le_x = 0.0
            try:
                le_x = float(airfoil_sdf.leading_edge[0])
            except AttributeError:
                pass
            le_wall = InsertWall(
                kind="leading_edge",
                skin_thickness=float(skin_thickness),
                x_limit=float(leading_edge_insert_x_m),
                le_x=le_x,
                le_y=0.0,
                radius=le_radius,
            )
            le = le_wall.to_sdf_component(
                airfoil_sdf, airfoil_expr_node=self._airfoil_expr
            )
            wall_recipes[_label_str(ComponentLabel.LEADING_EDGE_INSERT)] = le_wall
            self._components[ComponentLabel.LEADING_EDGE_INSERT] = le

        self._components_unrotated = dict(self._components)
        self._wall_recipes = wall_recipes
        self._skin_outer_boundary_unrotated_sdf = skin_outer_sdf
        self._skin_inner_boundary_unrotated_sdf = skin_inner_sdf
        if abs(self._twist) > 1e-10:
            skin_outer_sdf = rotate_field(skin_outer_sdf, self._twist)
            skin_inner_sdf = rotate_field(skin_inner_sdf, self._twist)
            twist_expr = Rotate(self._airfoil_expr, angle=self._twist, cx=0.0, cy=0.0)
            self._skin_outer_boundary_expr = Offset(twist_expr, amount=0.0)
            self._skin_inner_boundary_expr = Offset(twist_expr, amount=-skin_thickness)
            self._components = {
                lbl: twist_component(comp, self._twist)
                for lbl, comp in self._components_unrotated.items()
            }
        self._skin_outer_boundary_sdf = skin_outer_sdf
        self._skin_inner_boundary_sdf = skin_inner_sdf
        self._eval_union_expr = union_expr(
            [c.expr for c in self._components.values()]
        )

    def __getitem__(self, label):
        if label not in self._components:
            raise KeyError(f"Unknown component '{label}'. Available: {self.labels}")
        return self._components[label]

    def __iter__(self):
        return iter(self._components)

    def __len__(self):
        return len(self._components)

    @property
    def labels(self):
        return list(self._components.keys())

    @property
    def airfoil(self):
        return self._af

    @property
    def twist_rad(self) -> float:
        """Section twist angle in radians (S → B frame)."""
        return float(self._twist)

    @property
    def underlying_section(self):
        """Concrete section graph host (self for :class:`MultiCellSection`)."""
        return self

    @property
    def components(self) -> dict[str, Any]:
        """Active (possibly twisted) labelled component map."""
        return self._components

    @property
    def components_unrotated(self) -> dict[str, Any]:
        """Public chord-frame skin/web/core map (:class:`SectionMidlineHost`)."""
        return self._components_unrotated

    @property
    def spar_cap_components_unrotated(self) -> dict[str, Any]:
        """Public chord-frame spar-cap map (:class:`SectionMidlineHost`)."""
        return dict(self._spar_cap_components_unrotated)

    @property
    def wall_recipes(self) -> dict[str, Any]:
        """Frozen wall recipes keyed by label (assemble-time SoT)."""
        return dict(self._wall_recipes)

    @property
    def spar_cap_wall_recipes(self) -> dict[str, Any]:
        """Frozen discrete spar-cap recipes (chord frame; midline SoT)."""
        return dict(self._spar_cap_wall_recipes)

    @property
    def structural_family(self):
        """Catalogue family letter echo (not assemble control flow)."""
        return self._structural_family

    @property
    def cap_placement(self) -> CapPlacementRule:
        """Effective :class:`CapPlacementRule` used at assemble."""
        return self._cap_placement

    def eval(self, label, x, y, memo=None):
        return self._components[label](x, y, memo=memo)

    def eval_all(self, x, y):
        memo: dict = {}
        return {label: comp(x, y, memo=memo) for label, comp in self._components.items()}

    def eval_union(self, x, y):
        from blade_geometry_model.section_geometry.sdf.ir import eval_expr_batch

        return eval_expr_batch(self._eval_union_expr, x, y, {})

    @classmethod
    def d_spar(cls, airfoil_sdf, web_x=0.25, **kwargs):
        """Single-web section with caps (catalogue 1C / per_web_pads).

        Historical name; continuous-box (D) with one web is illegal (no 1D).
        Defaults to family ``C`` unless the caller overrides ``structural_family``
        or ``cap_placement``.
        """
        kwargs.setdefault("structural_family", "C")
        return cls(airfoil_sdf, web_x_positions=[web_x], **kwargs)

    @classmethod
    def twin_web(cls, airfoil_sdf, x_fore=0.20, x_aft=0.50, **kwargs):
        return cls(airfoil_sdf, web_x_positions=[x_fore, x_aft], **kwargs)

    @classmethod
    def torsion_box(cls, airfoil_sdf, x_fore=0.15, x_mid=0.35, x_aft=0.55, **kwargs):
        return cls(airfoil_sdf, web_x_positions=[x_fore, x_mid, x_aft], **kwargs)

    @property
    def n_webs(self):
        return sum(
            1
            for k in self._components
            if str(getattr(k, "value", k)).startswith(("shear_web_", "web_"))
        )

    @property
    def n_cells(self):
        return self.n_webs + 1

    def __repr__(self):
        return (
            f"MultiCellSection("
            f"cap_placement={self._cap_placement.kind!r}, "
            f"n_webs={self.n_webs}, "
            f"n_cells={self.n_cells}, "
            f"twist={np.degrees(self._twist):.1f}°, "
            f"components={self.labels})"
        )
