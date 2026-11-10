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

import warnings
from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any

from implicit_geometry_kernel.section.spec.runner_lifecycle import RunnerLifecycleMixin
from implicit_geometry_kernel.section.spec.options import SectionOptions
from implicit_geometry_kernel.section.spec.section_spec import SectionSpec
from implicit_geometry_kernel.section.spec.worker import SectionGeometryResult, StationLayoutResolution, build_section_geometry, resolve_station_layout

if TYPE_CHECKING:
    from implicit_geometry_kernel.catalogue.system_layout import SystemLayoutSpec
    from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF
    from implicit_geometry_kernel.section.layout.section_layout import SectionLayout


class SectionGeometry(RunnerLifecycleMixin):
    """Resolve the effective station layout, then build the section geometry.

    Parameters
    ----------
    spec:
        Frozen :class:`SectionSpec` (layout, airfoil, gauges, twist, webs).
    options:
        :class:`SectionOptions` (output gates, quadtree, sample counts).
    layout, airfoil:
        Deprecated keyword construction; use ``spec`` / ``options`` or
        :meth:`from_kwargs`. Nominal catalogue ``SystemLayoutSpec`` or
        parametric ``SectionLayout``, and station ``AirfoilSDF``. The
        *effective* layout may step down if the local cavity cannot anchor
        every web.
    skin_thickness, spar_cap_thickness, shear_web_thickness:
        Design-variable thicknesses [m].
    twist_rad:
        Geometric twist [rad] for the S→B map.
    pivot_fraction:
        Pitch-axis chord fraction from the LE (S-frame origin). Default: the
        standalone frame fallback ``0.333``.
    outputs:
        Compute-cost gates — ``properties`` / ``midlines`` / ``cell_wall_loops``.
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
        spec: SectionSpec | None = None,
        options: SectionOptions | None = None,
        *,
        layout: SectionLayout | SystemLayoutSpec | str | None = None,
        airfoil: AirfoilSDF | None = None,
        skin_thickness: float | None = None,
        spar_cap_thickness: float | None = None,
        shear_web_thickness: float | None = None,
        twist_rad: float | None = None,
        twist_angle_rad: float | None = None,
        twist_angle: float | None = None,
        outputs: Iterable[str] = ("properties", "midlines", "cell_wall_loops"),
        web_offsets: Sequence[float] | None = None,
        min_cavity_depth_m: float = 0.0,
        n_web_samples: int = 20,
        n_cap_samples: int = 80,
        retain_mesh: bool = False,
        props_backend: Any | None = None,
        quadtree_spec: Any | None = None,
        airfoil_for_spec: Any | None = None,
        resolution: StationLayoutResolution | None = None,
        spar_cap_ply_pitch_m: float | None = None,
        zone_gauges_m: Mapping[str, float] | None = None,
        pivot_fraction: float | None = None,
    ) -> None:
        from implicit_geometry_kernel.section.spec.section_spec import resolve_shear_web_thickness, resolve_skin_gauges, resolve_spar_cap_thickness, resolve_twist_rad

        if spec is not None:
            if layout is not None or airfoil is not None:
                raise TypeError("SectionGeometry() takes spec/options or layout kwargs, not both")
            self._bind_spec(
                spec,
                options or SectionOptions(),
                resolution=resolution,
                min_cavity_depth_m=min_cavity_depth_m,
            )
            return
        if layout is None or airfoil is None:
            raise TypeError(
                "SectionGeometry requires a SectionSpec or layout= and airfoil= kwargs"
            )
        warnings.warn(
            "SectionGeometry(layout=..., airfoil=...) is deprecated; "
            "use SectionGeometry(spec, options) or SectionGeometry.from_kwargs. "
            "Removed in the next minor.",
            DeprecationWarning,
            stacklevel=2,
        )
        self._layout = layout
        self._airfoil = airfoil
        self._zone_gauges_m = None if zone_gauges_m is None else dict(zone_gauges_m)
        if self._zone_gauges_m is not None:
            from implicit_geometry_kernel.section.spec.prep import _compat_kwargs_from_zone_gauges

            if any(
                v is not None
                for v in (skin_thickness, spar_cap_thickness, shear_web_thickness)
            ):
                raise ValueError(
                    "geometry received both zone_gauges_m and component thickness kwargs"
                )
            derived = _compat_kwargs_from_zone_gauges(self._zone_gauges_m)
            skin_thickness = derived["skin_thickness"]
            spar_cap_thickness = derived["spar_cap_thickness"]
            shear_web_thickness = derived["shear_web_thickness"]
        self._t_skin_gauges = resolve_skin_gauges(
            skin_thickness=skin_thickness, required=True
        )
        self._t_skin = max(self._t_skin_gauges)
        self._t_spar = resolve_spar_cap_thickness(
            spar_cap_thickness=spar_cap_thickness,
            default=0.0,
        )
        self._t_web = resolve_shear_web_thickness(
            shear_web_thickness=shear_web_thickness,
            default=0.0,
        )
        self._twist_rad = resolve_twist_rad(
            twist_rad=twist_rad,
            twist_angle_rad=twist_angle_rad,
            twist_angle=twist_angle,
        )
        self._outputs = tuple(outputs)
        self._web_x_positions_m = web_offsets
        self._min_cavity_depth_m = float(min_cavity_depth_m)
        self._n_web_samples = int(n_web_samples)
        self._n_cap_samples = int(n_cap_samples)
        self._retain_mesh = bool(retain_mesh)
        self._props_backend = props_backend
        self._quadtree_spec = quadtree_spec
        # Sizes the quadtree domain. Distinct from ``airfoil`` on purpose: the
        # QT box must be built on the B-frame (twisted) airfoil AABB, because a
        # chord-frame box clips rotated trailing-edge tips. Expressible here so
        # the runner is a full facade over ``build_section_geometry`` — a caller
        # that had to reach past it for one argument would just keep calling the
        # function.
        self._airfoil_for_spec = airfoil_for_spec
        self._executed = False
        # A caller that has already resolved this station passes the resolution
        # in rather than the three values derived from it (effective layout, web
        # x positions, web y coords). Those three can drift apart; a resolution
        # cannot disagree with itself.
        #
        # It also keeps a *frozen* topology frozen. A sizing loop resolves once
        # at the minimum skin so the design-variable set is fixed for the whole
        # search; re-resolving at each trial's thickness could hand back a
        # different web set than the station was frozen against.
        self._layout_frozen = resolution is not None
        self._resolution: StationLayoutResolution | None = resolution
        self._result: SectionGeometryResult | None = None
        self._spar_cap_ply_pitch_m = (
            float(spar_cap_ply_pitch_m)
            if spar_cap_ply_pitch_m is not None
            else None
        )
        self._pivot_fraction = (
            float(pivot_fraction) if pivot_fraction is not None else None
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
        self._t_spar = float(spec.gauges.cap)
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

    @classmethod
    def from_kwargs(cls, **kwargs: Any) -> SectionGeometry:
        """Deprecated: build from historical layout/airfoil keywords."""
        warnings.warn(
            "SectionGeometry.from_kwargs is deprecated; construct "
            "SectionGeometry(spec, options). Removed in the next minor.",
            DeprecationWarning,
            stacklevel=2,
        )
        layout = kwargs.pop("layout", None)
        airfoil = kwargs.pop("airfoil", None)
        if layout is None or airfoil is None:
            raise TypeError("from_kwargs requires layout= and airfoil=")
        resolution = kwargs.pop("resolution", None)
        min_cavity_depth_m = kwargs.pop("min_cavity_depth_m", 0.0)
        spec = SectionSpec.from_kwargs(layout, airfoil, **kwargs)
        options = SectionOptions.from_kwargs(**kwargs)
        return cls(
            spec,
            options,
            resolution=resolution,
            min_cavity_depth_m=min_cavity_depth_m,
        )

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
