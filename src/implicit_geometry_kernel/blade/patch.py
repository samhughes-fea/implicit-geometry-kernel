"""Smooth partition patches x(t, z) with analytic derivatives and curvature (decision 0017).

A :class:`PartitionPatch` rebuilds one partition of the shell (a skin partition
or a web) as a smooth surface from its run at a list of span stations:

* across the partition, each station's run is splined through its own points
  by normalised arc length ``t`` (not-a-knot cubic, C2), sampled at ``n_t``
  values of ``t``, and those samples are splined again in ``t``;
* along the span, each of those ``t`` rows is a PCHIP law in ``z`` (C1), checked
  for station density like any PCHIP law (decision 0014).

From the patch the kernel evaluates ``x``, ``x_t``, ``x_z``, ``x_tt``, ``x_tz``
and ``x_zz`` analytically (no finite differences), and from them the unit
normal, the principal curvatures ``k1 <= k2`` with their directions, and a
kind tag. :meth:`PartitionPatch.closest` projects points onto the patch by
Newton's method on (t, z), so the signed distance and its gradient (the normal)
near the surface come from the same analytic derivatives.

Accuracy is set by the run's own points: on a circular section the interior
curvature is within about 0.1 % of the mid-surface radius. Near a partition end,
a vertex the section inserted on a chord (a cap edge or web foot) can move the
end curvature by several per cent.

Sign convention: the normal is ``x_t x x_z`` normalised; a curvature is positive
when the surface bends towards the normal.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

from implicit_geometry_kernel.blade.skeleton import run_t
from implicit_geometry_kernel.structure.laws import PchipLaw, SamplingReport

KINDS = ("FLAT", "SINGLY", "DOUBLY", "UMBILIC")


def _spline_second(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Knot second derivatives of the not-a-knot cubic spline through ``(t, y)`` (C2).

    Not-a-knot ends keep the curvature at a partition's ends instead of forcing
    it to zero. ``y`` may be ``(n, m)``: each column is splined. Two knots give a
    line, three a parabola.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    n = len(t)
    m = np.zeros_like(y)
    if n < 3:
        return m
    h = np.diff(t)
    d = np.diff(y, axis=0) / (h[:, None] if y.ndim == 2 else h)
    if n == 3:
        m[:] = 2.0 * (d[1] - d[0]) / (h[0] + h[1])
        return m
    a = np.zeros((n, n))
    rhs = np.zeros_like(y)
    a[0, 0], a[0, 1], a[0, 2] = h[1], -(h[0] + h[1]), h[0]
    a[-1, -3], a[-1, -2], a[-1, -1] = h[-1], -(h[-2] + h[-1]), h[-2]
    for k in range(1, n - 1):
        a[k, k - 1] = h[k - 1]
        a[k, k] = 2.0 * (h[k - 1] + h[k])
        a[k, k + 1] = h[k]
        rhs[k] = 6.0 * (d[k] - d[k - 1])
    return np.linalg.solve(a, rhs)


def _spline_eval(t: np.ndarray, y: np.ndarray, m: np.ndarray, tq: float, order: int) -> np.ndarray:
    """Natural cubic spline (value, first or second derivative) at one ``tq``."""
    n = len(t)
    k = int(np.clip(np.searchsorted(t, tq, side="right") - 1, 0, n - 2))
    h = t[k + 1] - t[k]
    a = (t[k + 1] - tq) / h
    b = (tq - t[k]) / h
    if order == 0:
        return a * y[k] + b * y[k + 1] + ((a**3 - a) * m[k] + (b**3 - b) * m[k + 1]) * h * h / 6.0
    if order == 1:
        return (y[k + 1] - y[k]) / h + (-(3 * a * a - 1) * m[k] + (3 * b * b - 1) * m[k + 1]) * h / 6.0
    return a * m[k] + b * m[k + 1]


class _LinearRow:
    """Linear interpolation in z with its piecewise derivatives (ruled bays)."""

    def __init__(self, z: np.ndarray, values: np.ndarray) -> None:
        self.z, self.values = np.asarray(z, float), np.asarray(values, float)

    def at(self, z: float) -> float:
        return float(np.interp(z, self.z, self.values))

    def derivative(self, z, order: int = 1) -> np.ndarray:
        zz = np.asarray(z, float)
        if order == 2:
            return np.zeros(zz.shape)
        k = np.clip(np.searchsorted(self.z, zz, side="right") - 1, 0, self.z.size - 2)
        return (self.values[k + 1] - self.values[k]) / (self.z[k + 1] - self.z[k])


def _dedupe(pts: np.ndarray) -> np.ndarray:
    keep = np.concatenate([[True], np.linalg.norm(np.diff(pts, axis=0), axis=1) > 1e-12])
    return pts[keep]


@dataclass(frozen=True)
class SurfacePoint:
    """Analytic surface data at one (t, z).

    ``kind`` is ``FLAT`` (both curvatures within ``tol``), ``SINGLY`` (one),
    ``UMBILIC`` (equal and non-zero) or ``DOUBLY``. Curvatures within the
    tolerance are reported as exactly zero.
    """

    t: float
    z: float
    x: np.ndarray
    x_t: np.ndarray
    x_z: np.ndarray
    normal: np.ndarray
    k1: float
    k2: float
    d1: np.ndarray
    d2: np.ndarray
    theta: float
    kind: str


class PartitionPatch:
    """One partition (skin run or web) as a smooth surface over (t, z).

    Parameters
    ----------
    geo:
        A blade (``BladeGeometry``).
    name:
        The run name on ``section_at(z).shell_graph()`` (loop partitions and webs).
    zs:
        Increasing span stations where the run exists. PCHIP along z needs dense
        stations: each ``t`` row is checked and too-sparse stations raise
        ``UnderSampledSpanLaw`` unless ``sampling="flag"``.
    n_t:
        Samples across the partition (default 33).
    rtol:
        Station-density tolerance as a fraction of the patch's in-plane extent
        (the larger of its x and y ranges).
    curvature_tol:
        Curvatures below this magnitude (1/m) are exactly zero and set the kind.
    z_interp:
        ``"pchip"`` or ``"linear"`` along the span; default: the blade's own
        ``span_interp``. A linear patch is ruled between stations (zero
        spanwise curvature inside a bay) and is not density-checked.
    """

    def __init__(
        self,
        geo,
        name: str,
        zs: Sequence[float],
        *,
        n_t: int = 33,
        rtol: float = 5e-3,
        sampling: str = "raise",
        curvature_tol: float = 1e-6,
        z_interp: str | None = None,
    ) -> None:
        z = np.asarray([float(v) for v in zs])
        if z.ndim != 1 or z.size < 2 or np.any(np.diff(z) <= 0.0):
            raise ValueError("zs must be two or more increasing span stations")
        self.name = str(name)
        self.z = z
        self.t = np.linspace(0.0, 1.0, int(n_t))
        self.curvature_tol = float(curvature_tol)
        self.z_interp = str(z_interp or getattr(geo, "span_interp", "linear"))
        if self.z_interp not in ("linear", "pchip"):
            raise ValueError("z_interp must be 'linear' or 'pchip'")
        grid = np.empty((self.t.size, z.size, 3))
        for j, zj in enumerate(z):
            graph = geo.section_at(float(zj)).shell_graph()
            runs = {r.name: r for r in graph.loop + graph.webs}
            if self.name not in runs:
                raise ValueError(f"run {self.name!r} is not on the shell graph at z = {zj:g}")
            pts = _dedupe(np.asarray(runs[self.name].points, float)[:, :2])
            tt = run_t(pts)
            mm = _spline_second(tt, pts)
            for i, ti in enumerate(self.t):
                p = _spline_eval(tt, pts, mm, float(ti), 0)
                grid[i, j] = (p[0], p[1], float(zj))
        self.grid = grid
        self.sampling: list[SamplingReport] = []
        self._rows: list[tuple] = []
        extent = float(max(np.ptp(grid[..., 0]), np.ptp(grid[..., 1]), 1e-12))
        for i in range(self.t.size):
            laws: list[Any] = []
            for c in (0, 1):
                values = grid[i, :, c]
                if self.z_interp == "linear":
                    laws.append(_LinearRow(z, values))
                    continue
                law = PchipLaw(
                    z.tolist(), values.tolist(), name=f"{self.name} t={self.t[i]:.3f} {'xy'[c]}",
                    rtol=0.0, atol=rtol * extent, sampling=sampling,
                )
                self.sampling.append(law.sampling)
                laws.append(law)
            self._rows.append((laws[0], laws[1]))

    @property
    def sampling_ok(self) -> bool:
        """True when every t row passed the station-density check."""
        return all(r.ok for r in self.sampling)

    def _column(self, z: float, order: int) -> np.ndarray:
        """``(n_t, 3)`` z-derivative ``order`` of the t rows at ``z``."""
        out = np.empty((self.t.size, 3))
        for i, (lx, ly) in enumerate(self._rows):
            if order == 0:
                out[i] = (lx.at(z), ly.at(z), z)
            else:
                out[i] = (
                    float(lx.derivative([z], order)[0]),
                    float(ly.derivative([z], order)[0]),
                    1.0 if order == 1 else 0.0,
                )
        return out

    def derivatives(self, t: float, z: float) -> dict[str, np.ndarray]:
        """``x``, ``x_t``, ``x_z``, ``x_tt``, ``x_tz``, ``x_zz`` at (t, z), analytic."""
        tq = float(np.clip(t, 0.0, 1.0))
        zq = float(np.clip(z, self.z[0], self.z[-1]))
        c0, c1, c2 = (self._column(zq, k) for k in (0, 1, 2))
        m0, m1, m2 = (_spline_second(self.t, c) for c in (c0, c1, c2))
        return {
            "x": _spline_eval(self.t, c0, m0, tq, 0),
            "x_t": _spline_eval(self.t, c0, m0, tq, 1),
            "x_tt": _spline_eval(self.t, c0, m0, tq, 2),
            "x_z": _spline_eval(self.t, c1, m1, tq, 0),
            "x_tz": _spline_eval(self.t, c1, m1, tq, 1),
            "x_zz": _spline_eval(self.t, c2, m2, tq, 0),
        }

    def __call__(self, t: float, z: float) -> np.ndarray:
        """Point ``x(t, z)`` in the blade frame, metres."""
        return self.derivatives(t, z)["x"]

    def at(self, t: float, z: float) -> SurfacePoint:
        """Normal, principal curvatures and directions, and kind at (t, z)."""
        d = self.derivatives(t, z)
        xt, xz = d["x_t"], d["x_z"]
        n = np.cross(xt, xz)
        n = n / np.linalg.norm(n)
        e, f, g = xt @ xt, xt @ xz, xz @ xz
        l, m, nn = d["x_tt"] @ n, d["x_tz"] @ n, d["x_zz"] @ n
        first = np.array([[e, f], [f, g]])
        second = np.array([[l, m], [m, nn]])
        shape = np.linalg.solve(first, second)
        vals, vecs = np.linalg.eig(shape)
        vals, vecs = np.real(vals), np.real(vecs)
        order = np.argsort(vals)
        k = vals[order]
        uv = vecs[:, order]
        dirs = [uv[0, i] * xt + uv[1, i] * xz for i in range(2)]
        dirs = [v / np.linalg.norm(v) for v in dirs]
        tol = self.curvature_tol
        k = np.where(np.abs(k) < tol, 0.0, k)
        if k[0] == 0.0 and k[1] == 0.0:
            kind = "FLAT"
        elif k[0] == 0.0 or k[1] == 0.0:
            kind = "SINGLY"
        elif abs(k[0] - k[1]) < tol:
            kind = "UMBILIC"
        else:
            kind = "DOUBLY"
        tu = xt / np.linalg.norm(xt)
        theta = float(np.arctan2(np.dot(np.cross(tu, dirs[0]), n), np.dot(tu, dirs[0])))
        if kind in ("FLAT", "UMBILIC"):
            theta = 0.0
        return SurfacePoint(float(t), float(z), d["x"], xt, xz, n, float(k[0]), float(k[1]),
                            dirs[0], dirs[1], theta, kind)

    def closest(self, point: Sequence[float], *, t0: float = 0.5, z0: float | None = None,
                tol: float = 1e-12, max_iter: int = 50) -> tuple[float, float, float]:
        """Newton projection of ``point`` onto the patch: ``(t, z, signed distance)``.

        The distance is positive on the normal's side. The gradient of the
        distance at ``point`` is the unit normal at the foot.
        """
        p = np.asarray(point, float)
        t = float(t0)
        z = float(np.clip(p[2] if z0 is None else z0, self.z[0], self.z[-1]))
        for _ in range(max_iter):
            d = self.derivatives(t, z)
            r = d["x"] - p
            g = np.array([r @ d["x_t"], r @ d["x_z"]])
            h = np.array([
                [d["x_t"] @ d["x_t"] + r @ d["x_tt"], d["x_t"] @ d["x_z"] + r @ d["x_tz"]],
                [d["x_t"] @ d["x_z"] + r @ d["x_tz"], d["x_z"] @ d["x_z"] + r @ d["x_zz"]],
            ])
            step = np.linalg.solve(h, -g)
            t_new = float(np.clip(t + step[0], 0.0, 1.0))
            z_new = float(np.clip(z + step[1], self.z[0], self.z[-1]))
            done = abs(t_new - t) < tol and abs(z_new - z) < tol * max(1.0, abs(z))
            t, z = t_new, z_new
            if done:
                break
        sp = self.at(t, z)
        return t, z, float((p - sp.x) @ sp.normal)


__all__ = ["KINDS", "PartitionPatch", "SurfacePoint"]
