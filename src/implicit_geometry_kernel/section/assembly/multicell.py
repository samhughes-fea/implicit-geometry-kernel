"""
MultiCellSection: an N-web blade section with N+1 bays.

The section is assembled through ``structure``: :func:`~.plan.plan_section` resolves
the inputs, the builders in :mod:`.members` turn them into contributions (skin wall,
stacked caps, web branches, core and insert regions), and
:func:`~implicit_geometry_kernel.structure.assemble` combines them. The labelled
components the rest of the section pipeline reads are the builders' components, in
build order; twist is applied last.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from implicit_geometry_kernel.sdf.nodes import Offset, Rotate, union_expr
from implicit_geometry_kernel.sdf.kernels.transforms import rotate_field
from implicit_geometry_kernel.section.assembly import members as section_members
from implicit_geometry_kernel.section.assembly.plan import SectionMember, plan_section
from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
from implicit_geometry_kernel.section.assembly.labels import ComponentLabel
from implicit_geometry_kernel.section.walls.base import twist_component
from implicit_geometry_kernel.structure import BuilderRegistry, Structure, assemble


def section_registry() -> BuilderRegistry:
    """A registry of the default section builders (skin, webs, caps, cores, inserts)."""
    registry = BuilderRegistry()
    for builder in section_members.SECTION_BUILDERS:
        registry.register(builder())
    return registry


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
        registry: BuilderRegistry | None = None,
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

        ``registry`` (optional): the member builders to run; defaults to
        :func:`section_registry`. A new member type is one more builder.
        """
        plan = plan_section(
            airfoil_sdf, web_x_positions, web_thickness, web_alignment, cap_height,
            skin_thickness, twist_angle, web_y_coords, trailing_edge_insert_x_m,
            leading_edge_insert_x_m, le_radius, core_enabled, y_search, structural_family,
            fixed_cap_anchor, pitch_fraction_of_chord_from_le, fixed_cap_chord_half_width,
            discrete_cap_chord_half_width, cap_placement=cap_placement,
            web_endpoints=web_endpoints, stage_sink=stage_sink,
            spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
        )
        contributions = (registry or section_registry()).run(plan)
        self._structure = assemble(plan, contributions)
        self._plan = plan

        self._af = airfoil_sdf
        self._structural_family = plan.structural_family
        self._cap_placement = plan.cap_rule
        self._twist = plan.twist
        self._skin_t = plan.skin_thickness
        self._skin_t_ss = plan.t_ss
        self._skin_t_ps = plan.t_ps
        self._airfoil_expr = plan.airfoil_expr

        self._components: dict[Any, Any] = {}
        wall_recipes: dict[str, Any] = {}
        pads: dict[Any, Any] = {}
        cap_recipes: dict[str, Any] = {}
        for c in contributions:
            member = c.source
            if not isinstance(member, SectionMember):
                continue
            self._components[member.key] = member.component
            wall_recipes.update(member.recipes)
            if c.spec.kind == "cap":
                pads.update(member.extra)
                cap_recipes.update(member.recipes)
        self._spar_cap_components_unrotated = pads
        self._spar_cap_wall_recipes = cap_recipes
        self._components_unrotated = dict(self._components)
        self._wall_recipes = wall_recipes

        t_skin = plan.skin_thickness
        # Boundary callables come from the suction skin (full annulus at its gauge);
        # the boundary exprs use the thicker of the two skins.
        skin_ss = self._components_unrotated[ComponentLabel.SKIN_SUCTION_SURFACE]
        skin_outer_sdf = skin_ss.outer_sdf
        skin_inner_sdf = skin_ss.inner_sdf
        self._skin_outer_boundary_expr = Offset(plan.airfoil_expr, amount=0.0)
        self._skin_inner_boundary_expr = Offset(plan.airfoil_expr, amount=-t_skin)
        self._skin_outer_boundary_unrotated_sdf = skin_outer_sdf
        self._skin_inner_boundary_unrotated_sdf = skin_inner_sdf
        if abs(self._twist) > 1e-10:
            skin_outer_sdf = rotate_field(skin_outer_sdf, self._twist)
            skin_inner_sdf = rotate_field(skin_inner_sdf, self._twist)
            twist_expr = Rotate(plan.airfoil_expr, angle=self._twist, cx=0.0, cy=0.0)
            self._skin_outer_boundary_expr = Offset(twist_expr, amount=0.0)
            self._skin_inner_boundary_expr = Offset(twist_expr, amount=-t_skin)
            self._components = {
                lbl: twist_component(comp, self._twist)
                for lbl, comp in self._components_unrotated.items()
            }
        self._skin_outer_boundary_sdf = skin_outer_sdf
        self._skin_inner_boundary_sdf = skin_inner_sdf
        self._eval_union_expr = union_expr([c.expr for c in self._components.values()])

    @property
    def structure(self) -> Structure:
        """The assembled member model, in the chord frame (before twist)."""
        return self._structure

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

    def eval_union(self, x, y):
        from implicit_geometry_kernel.sdf.nodes import eval_expr_batch

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
