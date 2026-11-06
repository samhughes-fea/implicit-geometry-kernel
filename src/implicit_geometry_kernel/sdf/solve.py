"""In-plane Newton solves on ``SDFNode`` fields (the kernel as a geometry oracle)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from implicit_geometry_kernel.sdf.nodes.base import SDFNode


@dataclass(frozen=True)
class SolveResult:
    points: np.ndarray
    residual: np.ndarray
    converged: np.ndarray


def _cap(step: np.ndarray, max_step: float) -> np.ndarray:
    n = np.linalg.norm(step, axis=1, keepdims=True)
    scale = np.minimum(1.0, max_step / np.maximum(n, 1e-30))
    return step * scale


@dataclass(frozen=True)
class GeometrySolver:
    """Newton projection with a finite-difference gradient (``eps`` in metres)."""

    tolerance: float = 1e-6
    max_step: float = 0.02
    iterations: int = 12
    eps: float = 1e-6

    def _phi(self, field: SDFNode, pts: np.ndarray) -> np.ndarray:
        return field.evaluate(pts).phi

    def gradient(self, field: SDFNode, pts: np.ndarray, plane: bool = True) -> np.ndarray:
        n = pts.shape[0]
        dims = 2 if plane else 3
        offsets = np.eye(3)[:dims] * self.eps
        batch = np.vstack([pts + o for o in offsets] + [pts - o for o in offsets])
        v = self._phi(field, batch)
        g = np.zeros((n, 3))
        for d in range(dims):
            g[:, d] = (v[d * n:(d + 1) * n] - v[(dims + d) * n:(dims + d + 1) * n]) / (2 * self.eps)
        return g

    def project(self, field: SDFNode, pts: np.ndarray, plane: bool = True) -> SolveResult:
        p = np.array(pts, dtype=np.float64, copy=True)
        f = self._phi(field, p)
        for _ in range(self.iterations):
            todo = np.abs(f) > self.tolerance
            if not todo.any():
                break
            q = p[todo]
            g = self.gradient(field, q, plane)
            g2 = np.maximum(np.einsum("ij,ij->i", g, g), 1e-24)
            step = _cap(-(f[todo] / g2)[:, None] * g, self.max_step)
            trial = q + step
            ft = self._phi(field, trial)
            better = np.abs(ft) < np.abs(f[todo])
            idx = np.nonzero(todo)[0][better]
            p[idx] = trial[better]
            f[idx] = ft[better]
            if not better.any():
                break
        return SolveResult(p, np.abs(f), np.abs(f) <= self.tolerance)

    def intersect(self, a: SDFNode, b: SDFNode, pts: np.ndarray) -> SolveResult:
        """In-plane solve of ``phi_a = phi_b = 0``."""
        p = np.array(pts, dtype=np.float64, copy=True)
        fa, fb = self._phi(a, p), self._phi(b, p)
        res = np.maximum(np.abs(fa), np.abs(fb))
        for _ in range(self.iterations):
            todo = res > self.tolerance
            if not todo.any():
                break
            q = p[todo]
            ga = self.gradient(a, q)[:, :2]
            gb = self.gradient(b, q)[:, :2]
            det = ga[:, 0] * gb[:, 1] - ga[:, 1] * gb[:, 0]
            ok = np.abs(det) > 1e-12
            dx = np.zeros((q.shape[0], 3))
            dx[ok, 0] = (-fa[todo][ok] * gb[ok, 1] + fb[todo][ok] * ga[ok, 1]) / det[ok]
            dx[ok, 1] = (-fb[todo][ok] * ga[ok, 0] + fa[todo][ok] * gb[ok, 0]) / det[ok]
            dx = _cap(dx, self.max_step)
            trial = q + dx
            ta, tb = self._phi(a, trial), self._phi(b, trial)
            tres = np.maximum(np.abs(ta), np.abs(tb))
            better = ok & (tres < res[todo])
            idx = np.nonzero(todo)[0][better]
            p[idx] = trial[better]
            fa[idx], fb[idx], res[idx] = ta[better], tb[better], tres[better]
            if not better.any():
                break
        return SolveResult(p, res, res <= self.tolerance)

    def zero_along(self, field: SDFNode, pts: np.ndarray, dirs: np.ndarray,
                   slope: np.ndarray) -> SolveResult:
        p = np.asarray(pts, dtype=np.float64)
        d = np.asarray(dirs, dtype=np.float64)
        t = np.zeros(p.shape[0])
        s = np.asarray(slope, dtype=np.float64)
        for _ in range(self.iterations):
            phi = self._phi(field, p + t[:, None] * d)
            if np.max(np.abs(phi), initial=0.0) <= self.tolerance:
                break
            t -= phi / s
        out = p + t[:, None] * d
        res = np.abs(self._phi(field, out))
        return SolveResult(out, res, (res <= self.tolerance) & (t >= 0.0))
