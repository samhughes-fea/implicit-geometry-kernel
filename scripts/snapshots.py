"""Golden numbers and the public-API snapshot.

    python scripts/snapshots.py            # check: exit 1 and list what moved
    python scripts/snapshots.py --update   # rewrite tests/golden/*.json

The golden numbers cover every layout in the catalogue at three twists, plus a
lofted blade cut at five stations. Refactors must leave them unchanged; a PR that
changes geometry on purpose updates them and says why in the changelog.

The API snapshot lists every public name (top-level ``__all__``, ``contract`` and
``structure``)
with its call signature. Adding a name is fine; removing or changing one needs a
deprecation first (see CHANGELOG.md).
"""

from __future__ import annotations

import argparse
import inspect
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "tests" / "golden"
SECTIONS_FILE = GOLDEN / "sections.json"
BLADE_FILE = GOLDEN / "blade.json"
API_FILE = GOLDEN / "api.json"

TWISTS_DEG = (0.0, 15.0, 35.0)
GAUGES = {"skin_thickness": 0.004, "spar_cap_thickness": 0.012, "shear_web_thickness": 0.005}
PROPERTY_KEYS = ("area", "cx", "cy", "Ixx", "Iyy", "Ixy")
BLADE_CUTS = (0.0, 0.5, 2.0, 4.5, 6.0)

# Floating-point noise differs a little between platforms and library versions.
RTOL = 1e-9
ATOL = 1e-12


def _polyline_summary(points) -> dict:
    """Order-sensitive fingerprint of a polyline: size, ends, sums and length."""
    p = np.asarray(points, dtype=float)
    p = p.reshape(-1, 2) if p.ndim == 1 else p
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1) if len(p) > 1 else np.zeros(0)
    return {
        "n": int(len(p)),
        "first": p[0].tolist() if len(p) else [],
        "last": p[-1].tolist() if len(p) else [],
        "sum": p.sum(axis=0).tolist(),
        "sum_sq": (p**2).sum(axis=0).tolist(),
        "length": float(seg.sum()),
    }


def section_case(key: str, twist_deg: float) -> dict:
    from implicit_geometry_kernel import Airfoil
    from implicit_geometry_kernel.section import SectionGeometry, resolve_system_type

    airfoil = Airfoil.naca4("2412").to_sdf(chord=1.0)
    result = SectionGeometry(
        layout=resolve_system_type(key),
        airfoil=airfoil,
        twist_rad=math.radians(twist_deg),
        outputs=("properties", "midlines", "cell_wall_loops"),
        **GAUGES,
    ).execute().get_results()
    props = {
        label: {k: float(v[k]) for k in PROPERTY_KEYS}
        for label, v in sorted((result.properties_map or {}).items())
    }
    strips = [
        {"label": str(s.label), "kind": str(s.kind), **_polyline_summary(s.midline_twisted_pitch_axis_blade_frame)}
        for s in sorted(result.midline_strips or (), key=lambda s: (str(s.label), str(s.kind)))
    ]
    cells = [
        {"id": str(c.id), "area": float(c.area_m2), **_polyline_summary(c.polyline_twisted_pitch_axis_blade_frame)}
        for c in sorted(result.cell_wall_loops or (), key=lambda c: str(c.id))
    ]
    return {"properties": props, "strips": strips, "cells": cells}


def compute_sections() -> dict:
    from implicit_geometry_kernel import SYSTEM_TYPE_KEYS

    out = {}
    for key in SYSTEM_TYPE_KEYS:
        for twist in TWISTS_DEG:
            out[f"{key}@{twist:g}"] = section_case(key, twist)
    return out


def compute_blade() -> dict:
    """A tapered, twisted two-web blade cut at five stations."""
    from implicit_geometry_kernel import Airfoil, BladeGeometry, Station

    stations = [
        Station(0.0, 1.0, Airfoil.naca4("0024"), twist=12.0, pitch_axis=0.3),
        Station(3.0, 0.8, Airfoil.naca4("2418"), twist=5.0, pitch_axis=0.3),
        Station(6.0, 0.5, Airfoil.naca4("2412"), twist=-2.0, pitch_axis=0.3),
    ]
    blade = BladeGeometry(stations, "2D-F", skin_thickness=0.004, spar_cap_thickness=0.012,
                          shear_web_thickness=0.005)
    out = {}
    for z in BLADE_CUTS:
        cut = blade.section_at(z)
        out[f"z={z:g}"] = {
            "oml": _polyline_summary(blade.oml_polyline(z)),
            "walls": [
                {"label": str(w.label), "kind": str(w.kind), "thickness": float(w.thickness), **_polyline_summary(w.points)}
                for w in sorted(cut.walls, key=lambda w: (str(w.label), str(w.kind)))
            ],
            "cells": [
                {"area": float(c.area), **_polyline_summary(c.points)}
                for c in sorted(cut.cells, key=lambda c: float(c.area))
            ],
        }
    return out


def _signature(obj) -> str:
    try:
        return str(inspect.signature(obj))
    except (TypeError, ValueError):
        return ""


def compute_api() -> dict:
    import implicit_geometry_kernel as k
    import implicit_geometry_kernel.contract as contract
    import implicit_geometry_kernel.structure as structure

    out = {}
    contract_names = getattr(contract, "__all__", None) or [
        n for n in dir(contract) if not n.startswith("_") and inspect.isclass(getattr(contract, n))
    ]
    for module, names in ((k, k.__all__), (contract, contract_names), (structure, structure.__all__)):
        for name in sorted(names):
            obj = getattr(module, name)
            entry = {"kind": type(obj).__name__, "signature": _signature(obj) if callable(obj) else ""}
            if inspect.isclass(obj):
                entry["kind"] = "class"
                entry["members"] = {
                    m: _signature(v) for m, v in sorted(vars(obj).items())
                    if not m.startswith("_") and (callable(v) or isinstance(v, (property, classmethod, staticmethod)))
                }
            out[f"{module.__name__}.{name}"] = entry
    return out


def diff(expected, actual, path="") -> list[str]:
    """Human-readable differences, with numbers compared to RTOL/ATOL."""
    if isinstance(expected, dict) and isinstance(actual, dict):
        out = [f"{path}/{k}: missing" for k in expected.keys() - actual.keys()]
        out += [f"{path}/{k}: new" for k in actual.keys() - expected.keys()]
        for k in expected.keys() & actual.keys():
            out += diff(expected[k], actual[k], f"{path}/{k}")
        return out
    if isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            return [f"{path}: length {len(expected)} -> {len(actual)}"]
        return [d for i, (e, a) in enumerate(zip(expected, actual)) for d in diff(e, a, f"{path}[{i}]")]
    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)) and not isinstance(expected, bool):
        if math.isclose(expected, actual, rel_tol=RTOL, abs_tol=ATOL):
            return []
        return [f"{path}: {expected!r} -> {actual!r}"]
    return [] if expected == actual else [f"{path}: {expected!r} -> {actual!r}"]


def _write(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--update", action="store_true", help="rewrite the golden files")
    args = parser.parse_args()
    jobs = ((SECTIONS_FILE, compute_sections), (BLADE_FILE, compute_blade), (API_FILE, compute_api))
    if args.update:
        for path, fn in jobs:
            _write(path, fn())
            print(f"wrote {path.relative_to(ROOT)}")
        return 0
    failed = False
    for path, fn in jobs:
        moved = diff(json.loads(path.read_text()), fn())
        if path == API_FILE:  # new public names are allowed; removals and changes are not
            moved = [m for m in moved if not m.endswith(": new")]
        print(f"{path.relative_to(ROOT)}: {'unchanged' if not moved else f'{len(moved)} differences'}")
        for line in moved[:40]:
            print("  " + line)
        failed |= bool(moved)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
