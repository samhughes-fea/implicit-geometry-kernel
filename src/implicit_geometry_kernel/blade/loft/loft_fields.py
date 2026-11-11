"""SDF oracles, query and size fields for BladeGeometry."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from implicit_geometry_kernel.blade.loft import BladeGeometry

from implicit_geometry_kernel.blade.cut import _twist_offset
from implicit_geometry_kernel.blade.shell_graph import PartitionEdge


class LoftFieldsMixin:
    def oml_sdf(self: BladeGeometry, meta: int = 0):
        """Exact 2-D polygon SDF of ``oml_polyline(z)`` at each point's own ``z``."""
        from implicit_geometry_kernel.sdf.nodes.lofted import LoftedPolygonSDF

        return LoftedPolygonSDF(lambda z: self.oml_polyline(float(z))[:, :2], meta=meta)

    def cavity_sdf(self: BladeGeometry, meta: int = 0):
        """Inner mould line: ``oml + T``. ``T`` is a step (skin, plus cap inside footprints); no taper."""
        from implicit_geometry_kernel.blade.topology.fields import StackThickness
        from implicit_geometry_kernel.sdf.nodes.tagged import FieldOffsetNode

        return FieldOffsetNode(self.oml_sdf(meta), StackThickness(self), scale=1.0)

    def wall_mid_sdf(self: BladeGeometry, meta: int = 0):
        """Signed polygon SDF of the ``shell_graph`` skin loop at each point's ``z``.

        The graph loop is the spec (run interiors and partition vertices on
        ``phi = 0``). Webs are branches, not this contour. Colour is ``wall_meta``.
        """
        from implicit_geometry_kernel.contract import InfeasibleSectionGeometry
        from implicit_geometry_kernel.sdf.nodes.lofted import LoftedPolygonSDF

        def provider(z: float):
            try:
                return self.section_at(float(z)).shell_graph().loop_polygon()
            except InfeasibleSectionGeometry:
                return self.oml_polyline(float(z))[:, :2]

        return LoftedPolygonSDF(provider, meta=meta)

    def wall_meta(self: BladeGeometry, table=None):
        """Tagged wall-strip DAG: ``evaluate((x,y,z)).meta`` is the region at that ``z``."""
        from implicit_geometry_kernel.blade.topology.fields import CutWallMetaSDF
        from implicit_geometry_kernel.sdf.nodes.base import FeatureTable

        return CutWallMetaSDF(self, table if table is not None else FeatureTable())

    def thickness_law(self: BladeGeometry):
        """:class:`~implicit_geometry_kernel.structure.laws.ThicknessLaw` over ``(z, s)``."""
        from implicit_geometry_kernel.structure.laws import ThicknessLaw

        geo = self

        class _BladeThickness(ThicknessLaw):
            def __call__(self, z, s) -> np.ndarray:
                return geo.thickness_at(z, s)

        return _BladeThickness()

    def named_edges(self: BladeGeometry, z: float):
        """Partition edges and ply-drop ``s`` bounds at ``z`` (mesher feature curves)."""

        z = float(z)
        edges = list(self.section_at(z).shell_graph().partition_edges())
        chord, _, pitch, offset, (x, upper, lower) = self.interpolate(z)
        for drop in self.ply_drops:
            if not (float(drop.z0) - 1e-12 <= z <= float(drop.z1) + 1e-12):
                continue
            for s_frac, tag in ((drop.s0, "s0"), (drop.s1, "s1")):
                s_frac = float(np.clip(s_frac, 0.0, 1.0))
                i = int(np.clip(round(s_frac * (len(x) - 1)), 0, len(x) - 1))
                local = np.array([[x[i] - pitch * chord, 0.5 * (upper[i] + lower[i])]], float)
                xy = _twist_offset(local, self.interpolate(z)[1], offset)[0]
                pt = np.array([xy[0], xy[1], z], float)
                name = f"{drop.name}:{tag}"
                edges.append(PartitionEdge(z, name, name, pt))
        return tuple(edges)

    def size_fields(self: BladeGeometry, points) -> dict[str, np.ndarray]:
        """Distance, thickness and a polyline curvature estimate at 3-D points."""
        pts = np.asarray(points, float).reshape(-1, 3)
        out = self.query(pts)
        curv = np.zeros(len(pts))
        for i, p in enumerate(pts):
            poly = self.oml_polyline(float(p[2]))[:, :2]
            j = int(np.argmin(np.linalg.norm(poly - p[:2], axis=1)))
            a = poly[(j - 1) % len(poly)]
            b = poly[j]
            c = poly[(j + 1) % len(poly)]
            ab, bc = b - a, c - b
            lab = float(np.linalg.norm(ab))
            lbc = float(np.linalg.norm(bc))
            if lab < 1e-15 or lbc < 1e-15:
                continue
            cross = ab[0] * bc[1] - ab[1] * bc[0]
            curv[i] = abs(cross) / max(lab * lbc * 0.5 * (lab + lbc), 1e-15)
        out["curvature"] = curv
        return out

    def thickness_at(self: BladeGeometry, z, s) -> np.ndarray:
        """Skin thickness at span ``z`` and chord fraction ``s``, plus ply-drop ramps."""
        z = np.asarray(z, float)
        s = np.asarray(s, float)
        z, s = np.broadcast_arrays(z, s)
        out = np.empty(z.shape, dtype=float)
        for idx in np.ndindex(z.shape):
            t = float(self._skin(float(z[idx])))
            for drop in self.ply_drops:
                t += drop.extra(float(z[idx]), float(s[idx]))
            out[idx] = t
        return out

    def shell_graph(self: BladeGeometry, z: float, include_solids: bool = False):
        """Loftable shell graph at ``z``. Solids (cores/inserts) are opt-in."""
        return self.section_at(z).shell_graph(include_solids=include_solids)

    def key_points(self: BladeGeometry, z: float) -> dict[str, np.ndarray]:
        """Named points at ``z``: LE, TE, cap edges, web feet."""
        return self.section_at(z).key_points()

    def query(self: BladeGeometry, points) -> dict[str, np.ndarray]:
        """Batched thickness, owning member, local frame and OML distance at 3-D points."""
        pts = np.asarray(points, float).reshape(-1, 3)
        oml = self.oml_sdf().evaluate(pts)
        field = self.wall_meta()
        sample = field.evaluate(pts)
        table = field.table
        names, kinds = [], []
        for m in sample.meta.astype(int):
            if m <= 0:
                names.append("")
                kinds.append("")
            else:
                spec = table[int(m)]
                names.append(spec.name)
                kinds.append(spec.kind)
        s = []
        for p in pts:
            chord, _, pitch, _, _ = self.interpolate(float(p[2]))
            s.append(float(np.clip((p[0] + pitch * chord) / max(chord, 1e-15), 0.0, 1.0)))
        T = np.array([float(np.ravel(self.thickness_at(p[2], ss))[0]) for p, ss in zip(pts, s)])
        frames = []
        for p in pts:
            _, twist, _, _, _ = self.interpolate(float(p[2]))
            c, si = np.cos(np.radians(twist)), np.sin(np.radians(twist))
            frames.append(np.array([[c, -si, 0.0], [si, c, 0.0], [0.0, 0.0, 1.0]]))
        return {
            "thickness": T,
            "member": np.array(names, dtype=object),
            "kind": np.array(kinds, dtype=object),
            "distance": oml.phi,
            "frame": np.stack(frames, axis=0),
        }

    def write_oml(self: BladeGeometry, path, n_span: int = 41, n_chord: int = 81) -> str:
        """Watertight OML as STL or OBJ from a quad grid of span rings."""
        from implicit_geometry_kernel.blade.oml import write_oml as _write_oml

        return _write_oml(self, path, n_span=n_span, n_chord=n_chord)

    def write_mid_surfaces(
        self: BladeGeometry,
        path,
        n_span: int = 41,
        n_chord: int = 81,
        members=None,
        include_solids: bool = False,
    ) -> str:
        """One lofted patch per authored region / web, as STL or OBJ (groups)."""
        from implicit_geometry_kernel.blade.oml import write_mid_surfaces as _write_mids

        return _write_mids(
            self, path, n_span=n_span, n_chord=n_chord,
            members=members, include_solids=include_solids,
        )

    def validate(self: BladeGeometry) -> list[str]:
        """Named station/member checks: span order, thin walls, members outside the OML."""
        issues = []
        zs = [s.z for s in self.stations]
        if any(b <= a for a, b in zip(zs, zs[1:])):
            issues.append("stations must increase in z")
        for i, st in enumerate(self.stations):
            if st.chord <= 0:
                issues.append(f"station {i} at z={st.z:g}: non-positive chord")
        for member in self.transverse:
            if member.z < zs[0] or member.z > zs[-1]:
                issues.append(f"transverse {member.name!r} at z={member.z:g} leaves the blade")
            if member.thickness < 1e-6:
                issues.append(f"transverse {member.name!r}: too thin")
        if self.tip_closure_thickness() < 1e-6:
            issues.append(f"tip closure {self.tip_closure.name!r}: too thin")
        for stringer in self.stringers:
            if stringer.extent is not None:
                if stringer.extent[0] < zs[0] - 1e-9 or stringer.extent[1] > zs[-1] + 1e-9:
                    issues.append(f"stringer {stringer.name!r} extent leaves the blade")
            if stringer.thickness < 1e-6:
                issues.append(f"stringer {stringer.name!r}: too thin")
        for z in (zs[0], 0.5 * (zs[0] + zs[-1]), zs[-1]):
            try:
                cut = self.section_at(z)
            except Exception as exc:
                issues.append(f"station z={z:g}: {exc}")
                continue
            oml = self.oml_polyline(z)
            if not np.allclose(oml[0], oml[-1]):
                issues.append(f"station z={z:g}: OML is not closed")
            for wall in cut.walls:
                if wall.thickness < 1e-6:
                    issues.append(f"station z={z:g} member {wall.label}: too thin")
        return issues

    def snap_oml(self: BladeGeometry, points, solver=None):
        """Project points onto the OML in each point's own ``z``-plane."""
        from implicit_geometry_kernel.sdf.solve import GeometrySolver

        solver = solver or GeometrySolver(tolerance=1e-8, max_step=0.05)
        return solver.project(self.oml_sdf(), np.asarray(points, float), plane=True)

    def snap_wall_mid(self: BladeGeometry, points, solver=None):
        """Project onto the ``shell_graph`` loop field in each point's ``z``-plane."""
        from implicit_geometry_kernel.sdf.solve import GeometrySolver

        solver = solver or GeometrySolver(tolerance=1e-8, max_step=0.05)
        return solver.project(self.wall_mid_sdf(), np.asarray(points, float), plane=True)
