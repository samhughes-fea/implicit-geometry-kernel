"""Shear web subcomponent."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.section_geometry.sdf.csg import offset
from blade_geometry_model.section_geometry.sdf.ir import Capsule, Intersect, Offset

from blade_geometry_model.section_geometry.sections.clip_algebra import (
    airfoil_expr,
)
from blade_geometry_model.section_geometry.sections.labels import (
    ComponentLabel,
)
from .base import SectionComponent


class ShearWeb(SectionComponent):
    """Constant-thickness shear web as a shell around a capsule axis."""

    label = ComponentLabel.SHEAR_WEB

    def __init__(
        self,
        airfoil_sdf,
        skin_thickness,
        x_top,
        y_top,
        x_bot,
        y_bot,
        thickness,
        alignment="chord_normal",
        twist_angle=0.0,
        *,
        airfoil_expr_node=None,
    ):
        if alignment not in ("chord_normal", "flapwise"):
            raise ValueError("alignment must be 'chord_normal' or 'flapwise'.")

        self.x_top = float(x_top)
        self.y_top = float(y_top)
        self.x_bot = float(x_bot)
        self.y_bot = float(y_bot)
        self.thickness = float(thickness)
        self.alignment = alignment
        self.twist_angle = float(twist_angle)
        self.skin_thickness = float(skin_thickness)
        self._af = airfoil_sdf

        poly = airfoil_expr_node if airfoil_expr_node is not None else airfoil_expr(airfoil_sdf)
        if poly is None:
            raise ValueError("ShearWeb requires an AirfoilSDF with vertices for IR evaluation.")

        r = thickness / 2.0
        # Extend the capsule axis beyond the inner-skin anchors so the round
        # end-caps land outside the cavity. The intersection with the cavity
        # then yields a clean constant-width strip flush against the skin;
        # otherwise the disk ends bite asymmetrically into slanted skins and
        # the web body appears tilted (runtime-confirmed on twisted sections).
        dx = float(x_top) - float(x_bot)
        dy = float(y_top) - float(y_bot)
        seg_len = float(np.hypot(dx, dy))
        ext = 2.0 * float(thickness) + float(skin_thickness)
        if seg_len > 1e-12:
            ux, uy = dx / seg_len, dy / seg_len
        else:
            ux, uy = 0.0, 1.0
        self.x_top_ext = float(x_top) + ux * ext
        self.y_top_ext = float(y_top) + uy * ext
        self.x_bot_ext = float(x_bot) - ux * ext
        self.y_bot_ext = float(y_bot) - uy * ext
        expr = Intersect(
            Capsule(
                ax=self.x_top_ext,
                ay=self.y_top_ext,
                bx=self.x_bot_ext,
                by=self.y_bot_ext,
                r=float(r),
            ),
            Offset(poly, amount=-skin_thickness),
        )
        super().__init__(expr)

    def _extend_to_skin_mid_surface(self, anchor, direction):
        """Walk ``anchor`` outward along ``direction`` onto the skin mid-surface.

        The anchors are inner-mould-line intercepts — the right locus for
        building the *solid*, but half a skin gauge short of where the web's
        **mid-surface** actually ends. Prefers the shared offset-ring oracle
        (:func:`midlines.junction_primitives.solve_web_foot_on_skin_mid`) so
        production feet and the topology oracle use the same polygon join
        contract; falls back to SDF ray-march on the mid-surface field.

        Never returns an unverified blind offset. A foot that cannot be placed
        on the mid-surface raises :class:`InfeasibleSectionGeometry` (SG-002)
        instead of handing ``cell_wall_loops`` a millimetre-scale miss.
        """
        from blade_geometry_model.section_contract.errors import (
            InfeasibleSectionGeometry,
        )

        half_gauge = 0.5 * float(self.skin_thickness)
        if half_gauge <= 0.0:
            return np.asarray(anchor, dtype=float)

        step = np.asarray(direction, dtype=float)
        norm = float(np.linalg.norm(step))
        if norm <= 1e-15:
            return np.asarray(anchor, dtype=float)
        step = step / norm
        origin = np.asarray(anchor, dtype=float)
        # Span well past IML and OML so tip / twist / stale-anchor cases still
        # register as interior segment hits on the mid ring.
        search = max(8.0 * half_gauge, 4.0 * float(self.skin_thickness), 1e-3)

        def _pick_outward_hit(hits):
            scored = []
            for h in hits:
                hh = np.asarray(h, dtype=float).reshape(2)
                along = float(np.dot(hh - origin, step))
                if along < -1e-12:
                    continue
                scored.append((along, float(np.linalg.norm(hh - origin)), hh))
            if not scored:
                return None
            scored.sort(key=lambda t: (t[0], t[1]))
            return scored[0][2]

        try:
            from blade_geometry_model.section_geometry.midlines.junction_primitives import (
                solve_web_foot_on_skin_mid,
            )

            axis_p0 = origin + search * step
            axis_p1 = origin - search * step
            hits = solve_web_foot_on_skin_mid(
                axis_p0, axis_p1, self._af, float(self.skin_thickness)
            )
            picked = _pick_outward_hit(hits) if hits else None
            if picked is not None:
                return picked
        except Exception:
            pass

        phi_mid = offset(self._af, -half_gauge)

        def phi_at(t):
            p = origin + float(t) * step
            return float(phi_mid(np.array([p[0]]), np.array([p[1]]))[0])

        # Sample the axis both ways; any sign-change bracket is a valid foot.
        samples = np.linspace(-search, search, 65)
        phis = np.asarray([phi_at(t) for t in samples], dtype=float)
        bracket = None
        for i in range(len(samples) - 1):
            p0, p1 = float(phis[i]), float(phis[i + 1])
            if p0 == 0.0:
                return origin + float(samples[i]) * step
            if p0 * p1 > 0.0:
                continue
            # Prefer the outward (t >= 0) crossing nearest the origin.
            t_lo, t_hi = float(samples[i]), float(samples[i + 1])
            if t_hi < 0.0:
                continue
            lo, hi = t_lo, t_hi
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if phi_at(mid) <= 0.0:
                    lo = mid
                else:
                    hi = mid
            t_hit = 0.5 * (lo + hi)
            if t_hit < -1e-12:
                continue
            if bracket is None or t_hit < bracket[0]:
                bracket = (t_hit, origin + t_hit * step)
        if bracket is not None:
            return bracket[1]

        raise InfeasibleSectionGeometry(
            "ShearWeb foot cannot reach the skin mid-surface along its axis "
            f"(anchor={origin.tolist()}, skin_thickness={float(self.skin_thickness):.6g}). "
            "Frozen web_y anchors are usually stale relative to the current skin, "
            "or the local cavity is thinner than the skin gauge (SG-002)."
        )

    def midline_polyline(self, n=20):
        n_use = max(2, int(n))
        top = np.array([self.x_top, self.y_top], dtype=float)
        bot = np.array([self.x_bot, self.y_bot], dtype=float)
        axis = top - bot
        return np.linspace(
            self._extend_to_skin_mid_surface(top, axis),
            self._extend_to_skin_mid_surface(bot, -axis),
            n_use,
        )

    def __repr__(self):
        return (
            f"ShearWeb(x=[{self.x_top:.3g},{self.x_bot:.3g}], "
            f"alignment={self.alignment!r}, "
            f"twist={np.degrees(self.twist_angle):.2f}°, "
            f"thickness={self.thickness:.3g})"
        )
