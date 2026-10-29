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

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from .worker import (
    SectionGeometryResult,
    StationLayoutResolution,
    build_section_geometry,
    resolve_station_layout,
    section_geometry_prep_key,
)


class SectionGeometryRunner:
    """Resolve the effective station layout, then build the section geometry.

    Parameters
    ----------
    layout:
        Nominal catalogue ``SystemLayoutSpec`` or parametric ``SectionLayout``
        requested for this station. The *effective* layout may step down if the
        local cavity cannot anchor every web (catalogue rung, or keyless
        ``N`` webs + ``CapPlacementRule``).
    airfoil:
        Station ``AirfoilSDF``.
    skin_thickness, spar_cap_thickness, shear_web_thickness:
        Design-variable thicknesses [m].
    twist_rad:
        Geometric twist [rad] for the S→B map.
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
        *,
        layout: Any,
        airfoil: Any,
        skin_thickness: float | None = None,
        spar_cap_thickness: float | None = None,
        shear_web_thickness: float | None = None,
        twist_rad: float | None = None,
        twist_angle_rad: float | None = None,
        twist_angle: float | None = None,
        outputs: Iterable[str] = ("properties", "midlines", "cell_wall_loops"),
        web_x_positions_m: Sequence[float] | None = None,
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
    ) -> None:
        from .kwarg_aliases import (
            resolve_shear_web_thickness,
            resolve_skin_gauges,
            resolve_spar_cap_thickness,
            resolve_twist_rad,
        )

        if layout is None:
            raise ValueError("layout is required")
        if airfoil is None:
            raise ValueError("airfoil is required")
        self._layout = layout
        self._airfoil = airfoil
        self._zone_gauges_m = None if zone_gauges_m is None else dict(zone_gauges_m)
        if self._zone_gauges_m is not None:
            from .worker import _compat_kwargs_from_zone_gauges

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
        self._web_x_positions_m = web_x_positions_m
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
        self._resolution: StationLayoutResolution | None = resolution
        self._result: SectionGeometryResult | None = None
        self._spar_cap_ply_pitch_m = (
            float(spar_cap_ply_pitch_m)
            if spar_cap_ply_pitch_m is not None
            else None
        )

    def resolve(self) -> "SectionGeometryRunner":
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

    def execute(self) -> "SectionGeometryRunner":
        """Resolve layout then build geometry. Idempotent; returns ``self``."""
        if self._executed:
            return self
        self.resolve()
        resolution = self._resolution
        assert resolution is not None  # set by resolve()
        # Build against the *effective* layout: the nominal one may name webs
        # this station cannot anchor, and sizing must not allocate design
        # variables for components the section does not have.
        self._result = build_section_geometry(
            resolution.layout,
            self._airfoil,
            **(
                {"zone_gauges_m": self._zone_gauges_m}
                if self._zone_gauges_m is not None
                else {
                    "skin_thickness": self._t_skin_gauges,
                    "spar_cap_thickness": self._t_spar,
                    "shear_web_thickness": self._t_web,
                }
            ),
            twist_rad=self._twist_rad,
            outputs=self._outputs,
            props_backend=self._props_backend,
            n_web_samples=self._n_web_samples,
            n_cap_samples=self._n_cap_samples,
            retain_mesh=self._retain_mesh,
            quadtree_spec=self._quadtree_spec,
            airfoil_for_spec=self._airfoil_for_spec,
            web_x_positions_m=resolution.web_x_positions_m,
            web_y_coords=getattr(resolution, "web_y_coords", None),
            spar_cap_ply_pitch_m=self._spar_cap_ply_pitch_m,
        )
        self._resolution = resolution
        self._executed = True
        return self

    def get_results(self) -> SectionGeometryResult:
        """Return the built :class:`SectionGeometryResult`. Requires :meth:`execute`."""
        if not self._executed or self._result is None:
            raise RuntimeError(
                "SectionGeometryRunner.execute() must be called before get_results()"
            )
        return self._result

    def get_layout_resolution(self) -> StationLayoutResolution:
        """Return the effective layout (anchored webs, drops).

        Requires :meth:`resolve` **or** :meth:`execute` — the resolution is the
        cheap first product and does not depend on the build.
        """
        if self._resolution is None:
            raise RuntimeError(
                "SectionGeometryRunner.resolve() or .execute() must be called "
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
            shear_web_thickness=self._t_web,
            twist_rad=self._twist_rad,
            outputs=self._outputs,
            web_x_positions_m=web_x,
            zone_gauges_m=self._zone_gauges_m,
        )


__all__ = ["SectionGeometryRunner"]
