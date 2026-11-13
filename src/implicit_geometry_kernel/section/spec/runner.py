"""Geometry worker runner: layout resolution + section build.

One lifecycle seam for the geometry stage, per
`docs/conventions/api/API_STANDARDS.md` § Worker runners. Wraps
:func:`resolve_station_layout` and :func:`build_section_geometry`; both stay
exported and remain the unit-test surface.

The two are sequenced deliberately: the layout a station can actually carry is
resolved *before* anything is assembled, so a dropped web is known to the caller
without paying for a build. ``get_layout_resolution()`` exposes that first
product; ``get_results()`` returns the built geometry.
"""

from __future__ import annotations

from implicit_geometry_kernel.section.spec.runner_lifecycle import RunnerLifecycleMixin
from implicit_geometry_kernel.section.spec.options import SectionOptions
from implicit_geometry_kernel.section.spec.section_spec import SectionSpec
from implicit_geometry_kernel.section.spec.worker import SectionGeometryResult, StationLayoutResolution, build_section_geometry, resolve_station_layout


class SectionGeometry(RunnerLifecycleMixin):
    """Resolve the effective station layout, then build the section geometry.

    Parameters
    ----------
    spec:
        Frozen :class:`SectionSpec` (layout, airfoil, gauges, twist, webs).
        The *effective* layout may step down if the local cavity cannot
        anchor every web.
    options:
        :class:`SectionOptions` (output gates, quadtree, sample counts).
    min_cavity_depth_m:
        Forwarded to layout resolution: a web whose local cavity is too shallow
        to build minimum-gauge caps against is dropped rather than anchored.

    Raises
    ------
    InfeasibleSectionGeometry
        From layout resolution when no web can be anchored — an over-constrained
        station, which callers must distinguish from a merely too-thick trial.
    """

    def __init__(
        self,
        spec: SectionSpec,
        options: SectionOptions | None = None,
        *,
        min_cavity_depth_m: float = 0.0,
        resolution: StationLayoutResolution | None = None,
    ) -> None:
        self._bind_spec(
            spec,
            options or SectionOptions(),
            resolution=resolution,
            min_cavity_depth_m=min_cavity_depth_m,
        )

    def _bind_spec(
        self,
        spec: SectionSpec,
        options: SectionOptions,
        *,
        resolution: StationLayoutResolution | None = None,
        min_cavity_depth_m: float = 0.0,
    ) -> None:
        self._spec = spec
        self._options = options
        self._layout = spec.layout
        self._airfoil = spec.airfoil
        self._zone_gauges_m = None if spec.zone_gauges_m is None else dict(spec.zone_gauges_m)
        self._t_skin_gauges = (spec.gauges.skin_ss, spec.gauges.skin_ps)
        self._t_skin = spec.gauges.skin_max
        cap = spec.gauges.cap
        self._t_spar = cap if isinstance(cap, tuple) else float(cap)
        self._t_web = spec.gauges.web_arg
        self._twist_rad = spec.twist_rad
        self._outputs = tuple(options.outputs)
        self._web_x_positions_m = spec.webs
        self._min_cavity_depth_m = float(min_cavity_depth_m)
        self._n_web_samples = int(options.n_web_samples)
        self._n_cap_samples = int(options.n_cap_samples)
        self._retain_mesh = bool(options.retain_mesh)
        self._props_backend = options.props_backend
        self._quadtree_spec = options.quadtree_spec
        self._airfoil_for_spec = options.airfoil_for_spec
        self._executed = False
        self._layout_frozen = resolution is not None
        self._resolution: StationLayoutResolution | None = resolution
        self._result: SectionGeometryResult | None = None
        self._spar_cap_ply_pitch_m = spec.spar_cap_ply_pitch_m
        self._pivot_fraction = spec.pivot_fraction

    def resolve(self) -> "SectionGeometry":
        """Resolve the effective layout **without building geometry**.

        The cheap half of the stage: which webs this station can actually anchor,
        and which were dropped. Idempotent; returns ``self``.

        Exists because callers legitimately want the layout far earlier than the
        geometry — a sizing loop resolves it once per station to fix the design
        variable set, then builds geometry many times per resolution while
        searching thicknesses. On the reference case that is 52 resolutions
        against 168 builds, so folding the two together would add roughly a
        third more builds for nothing.

        Until this existed the module docstring's promise that a dropped web is
        known "without paying for a build" was not true of any code path:
        ``execute()`` did both and ``get_layout_resolution()`` required it.
        """
        if self._resolution is None:
            self._resolution = resolve_station_layout(
                self._layout,
                self._airfoil,
                skin_thickness=self._t_skin,
                web_x_positions_m=self._web_x_positions_m,
                min_cavity_depth_m=self._min_cavity_depth_m,
            )
        return self

    def _assemble(self, resolution: StationLayoutResolution, *, live_y: bool = False):
        from implicit_geometry_kernel.section.spec import SectionOptions, SectionSpec

        y = None if live_y else getattr(resolution, "web_y_coords", None)
        thickness = (
            {"zone_gauges_m": self._zone_gauges_m}
            if self._zone_gauges_m is not None
            else {
                "skin_thickness": self._t_skin_gauges,
                "spar_cap_thickness": self._t_spar,
                "shear_web_thickness": self._t_web,
            }
        )
        spec = SectionSpec.from_kwargs(
            resolution.layout,
            self._airfoil,
            **thickness,
            twist_rad=self._twist_rad,
            web_offsets=resolution.web_x_positions_m,
            web_indices=resolution.web_indices or None,
            web_y_coords=y,
            spar_cap_ply_pitch_m=self._spar_cap_ply_pitch_m,
            pivot_fraction=self._pivot_fraction,
        )
        opts = SectionOptions.from_kwargs(
            outputs=self._outputs,
            props_backend=self._props_backend,
            n_web_samples=self._n_web_samples,
            n_cap_samples=self._n_cap_samples,
            retain_mesh=self._retain_mesh,
            quadtree_spec=self._quadtree_spec,
            airfoil_for_spec=self._airfoil_for_spec,
        )
        return build_section_geometry(spec, options=opts)


__all__ = ["SectionGeometry"]
