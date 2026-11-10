"""Execute / results lifecycle for SectionGeometry."""

from __future__ import annotations

from typing import Any, TYPE_CHECKING

from implicit_geometry_kernel.section.spec.layout_resolution import StationLayoutResolution
from implicit_geometry_kernel.section.spec.worker import SectionGeometryResult, section_geometry_prep_key

if TYPE_CHECKING:
    from implicit_geometry_kernel.section.spec.runner import SectionGeometry


class RunnerLifecycleMixin:
    _resolution: StationLayoutResolution | None
    _executed: bool
    _result: SectionGeometryResult | None

    def _assemble_dropping_webs(self, resolution: StationLayoutResolution, last_exc):
        """Wave 4: step down webs when flapwise re-anchor fails (not a frozen sizing freeze)."""
        from implicit_geometry_kernel.contract import InfeasibleSectionGeometry
        from implicit_geometry_kernel.section.spec.layout_resolution import (
            StationLayoutResolution,
            _reattach_catalogue_preset,
            _step_down_section_layout,
        )
        from implicit_geometry_kernel.section.spec.results import DroppedWeb

        xs = list(resolution.web_x_positions_m)
        idx = list(resolution.web_indices or range(len(xs)))
        dropped = list(resolution.dropped)
        while len(xs) > 1:
            fail_i = getattr(last_exc, "web_index", None)
            fail_x = getattr(last_exc, "web_x", None)
            if fail_i is not None and fail_i in idx:
                pos = idx.index(fail_i)
            elif fail_x is not None:
                pos = min(range(len(xs)), key=lambda j: abs(xs[j] - float(fail_x)))
            else:
                pos = 0
            lost = float(xs.pop(pos))
            lost_i = int(idx.pop(pos))
            dropped.append(
                DroppedWeb(
                    lost,
                    "flapwise web cannot re-anchor under twist; dropped at this station",
                    index=lost_i,
                )
            )
            sl = _step_down_section_layout(
                resolution.nominal_layout, len(xs), indices=idx,
            )
            sl = _reattach_catalogue_preset(sl, resolution.nominal_layout)
            new_res = StationLayoutResolution(
                sl, tuple(xs), tuple(dropped), resolution.nominal_layout,
                web_y_coords=(), web_indices=tuple(idx),
            )
            try:
                result = self._assemble(new_res, live_y=True)
                self._resolution = new_res
                return result
            except InfeasibleSectionGeometry as exc:
                last_exc = exc
        raise last_exc

    def execute(self) -> "SectionGeometry":
        """Resolve layout then build geometry. Idempotent; returns ``self``."""
        from implicit_geometry_kernel.contract import InfeasibleSectionGeometry

        if self._executed:
            return self
        self.resolve()
        resolution = self._resolution
        assert resolution is not None  # set by resolve()
        # Build against the *effective* layout: the nominal one may name webs
        # this station cannot anchor, and sizing must not allocate design
        # variables for components the section does not have.
        try:
            self._result = self._assemble(resolution)
        except InfeasibleSectionGeometry as exc:
            if self._layout_frozen or "re-anchor" not in str(exc).lower():
                raise
            try:
                self._result = self._assemble(resolution, live_y=True)
            except InfeasibleSectionGeometry:
                self._result = self._assemble_dropping_webs(resolution, exc)
        self._executed = True
        return self

    def get_results(self) -> SectionGeometryResult:
        """Return the built :class:`SectionGeometryResult`. Requires :meth:`execute`."""
        if not self._executed or self._result is None:
            raise RuntimeError(
                "SectionGeometry.execute() must be called before get_results()"
            )
        return self._result

    def get_layout_resolution(self) -> StationLayoutResolution:
        """Return the effective layout (anchored webs, drops).

        Requires :meth:`resolve` **or** :meth:`execute` — the resolution is the
        cheap first product and does not depend on the build.
        """
        if self._resolution is None:
            raise RuntimeError(
                "SectionGeometry.resolve() or .execute() must be called "
                "before get_layout_resolution()"
            )
        return self._resolution

    @property
    def prep_key(self) -> tuple[Any, ...]:
        """Cache key for this station's geometry inputs.

        Uses the **anchored** web positions once executed, falling back to the
        requested ones before that. Anchored positions change the built section
        without changing the layout key or any thickness, so a key built from the
        requested positions alone would let two different sections collide.
        """
        web_x = (
            self._resolution.web_x_positions_m
            if self._resolution is not None
            else self._web_x_positions_m
        )
        return section_geometry_prep_key(
            self._layout,
            self._airfoil,
            skin_thickness=self._t_skin,
            spar_cap_thickness=self._t_spar,
            shear_web_thickness=self._t_web,  # type: ignore
            twist_rad=self._twist_rad,
            outputs=self._outputs,
            web_x_positions_m=web_x,
            zone_gauges_m=self._zone_gauges_m,
        )
