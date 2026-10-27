"""A1 — what does the per-station chord-fraction web model actually do spanwise?

Web chordwise positions are resolved per station as ``frac x local_chord``
(``api/factory.py``), so every web migrates inboard as the chord tapers. Any web
that then fails to find inner-skin intercepts is walked *further* LE-ward by
``resolve_web_anchors`` -- up to 12 steps -- before being silently dropped.

This probe measures that, station by station, on the real reference rotor
geometry, for one layout per family that carries spar caps (B, C, D). It changes
nothing; it is the input to the root-anchored redesign.

Reported per station and per web:

* requested x (``frac x local_chord``) vs resolved anchor x;
* how far the web was walked LE-ward, absolute and as a fraction of chord;
* whether it was dropped, and how many survived;
* the spar-cap band the surviving webs imply, and the height
  ``_cap_feasibility_clamp`` allows -- which is family-specific:
  D keys on the depth at the aft end of the ``xs[0]..xs[-1]`` band,
  C clamps one discrete band per web, B is web-independent.

From the repo root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/\\
        spanwise_web_anchors.py
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

from _ensure_repo_path import ensure_repo_path  # noqa: E402  (script-local helper)

_REPO = ensure_repo_path()

from spanwise_web_anchors_shared import build_stations as _build_spanwise_stations  # noqa: E402
from blade_structure_precompute.module.section_selection.processing.adapters.scale_adapter import (  # noqa: E402
    airfoil_sdf_from_station,
)
from blade_geometry_model.section_geometry import (  # noqa: E402
    resolve_system_type,
)
from blade_geometry_model.section_geometry.sections.web_layout import (  # noqa: E402
    inner_y_at_x,
    resolve_web_anchors,
)
from blade_geometry_model.section_geometry.sections.spar_family import (  # noqa: E402
    _cap_feasibility_clamp,
)

logging.disable(logging.INFO)

ROTOR = "new_3.2m_rotor"
MATERIALS = {"skin": 17, "spar": 17, "web": 17}
# One cap-bearing layout per family. Family A has no spar caps, so the cap half
# of this probe does not apply there.
LAYOUTS = ("5D-F", "5C-CN", "5B-F")
# The reference sizing thicknesses; the point is the geometry, not the sizing.
T_SKIN = 0.001
T_SPAR = 0.004
T_WEB = 0.0015
OUT_DIR = Path(__file__).resolve().parents[1] / "output" / "spanwise_web_anchors"


def build_stations() -> list[dict[str, Any]]:
    """Return spanwise station data via the shared loader helper."""
    return _build_spanwise_stations(
        repo_root=_REPO,
        rotor=ROTOR,
        materials=MATERIALS,
        structural_system_key="5D-F",
    )


def probe_layout(  # noqa: C901 - a measurement script, linear by design
    layout_key: str, stations: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    layout = resolve_system_type(layout_key)
    family = str(layout.structural_family)
    fracs = tuple(float(f) for f in layout.web_chord_fracs)

    out: list[dict[str, Any]] = []
    for st in stations:
        rec: dict[str, Any] = {
            "idx": st["idx"],
            "zL": st["zL"],
            "chord": st["chord"],
            "n_webs_nominal": len(fracs),
        }
        try:
            airfoil = airfoil_sdf_from_station(st["coords"], st["metadata"])
        except Exception as exc:  # geometry could not be built at all
            rec["error"] = f"{type(exc).__name__}: {exc}"
            out.append(rec)
            continue

        # Exactly what api/factory.py does: frac * the *local* chord.
        xs_requested = [f * float(airfoil.chord) for f in fracs]
        try:
            anchors = resolve_web_anchors(airfoil, T_SKIN, xs_requested)
        except Exception as exc:
            rec["error"] = f"resolve_web_anchors: {type(exc).__name__}: {exc}"
            out.append(rec)
            continue

        xs_resolved = [float(a[0]) for a in anchors]
        rec["xs_requested"] = xs_requested
        rec["xs_resolved"] = xs_resolved
        rec["n_webs_resolved"] = len(xs_resolved)
        rec["n_dropped"] = len(fracs) - len(xs_resolved)

        # Pair each survivor with its nearest request to estimate the LE-ward walk.
        walks = []
        for xr in xs_resolved:
            nearest = min(xs_requested, key=lambda xq: abs(xq - xr))
            walks.append(float(nearest - xr))  # positive = walked LE-ward
        rec["walk_m"] = walks
        rec["walk_frac_chord"] = [w / max(float(airfoil.chord), 1e-12) for w in walks]
        rec["max_walk_frac_chord"] = max(walks, default=0.0) / max(float(airfoil.chord), 1e-12)

        # Cap band + allowance, per family coupling.
        rec["cap"] = cap_probe(family, airfoil, xs_resolved, layout)

        # Counterfactual: what the root-anchored model would do here. Positions
        # are the ROOT station's absolute offsets, reused unchanged. A web
        # survives only if the local cavity can actually hold it *at that x* --
        # the same inner-skin intercept test resolve_web_anchors uses, but with
        # the LE-ward walk removed, which is the whole point of the new model.
        #
        # Note "x <= chord" is NOT the right criterion: a web at 91% chord is
        # inside the section but in a trailing region too thin to anchor.
        xs_root = [f * float(stations[0]["chord"]) for f in fracs]
        survivors, dropped = [], []
        for x in xs_root:
            try:
                inner_y_at_x(airfoil, T_SKIN, float(x))
            except Exception as exc:
                dropped.append({"x_m": float(x), "why": type(exc).__name__})
                continue
            survivors.append(float(x))
        rec["root_anchored"] = {
            "xs_root_m": xs_root,
            "n_survivors": len(survivors),
            "n_dropped": len(dropped),
            "xs_survivors_m": survivors,
            "dropped": dropped,
        }
        if survivors and bool(layout.has_spar_caps):
            rec["root_anchored"]["cap"] = cap_probe(family, airfoil, survivors, layout)

        # How aggressive must dropping be to get a *usable* cap? Keeping the
        # thinnest-point sizing rule, drop the aft-most web repeatedly until the
        # allowance clears the material minimum, or until the family floor stops
        # us. This is the criterion A3 needs, calibrated rather than guessed.
        rec["greedy"] = greedy_drop_to_usable(family, airfoil, survivors, layout)
        out.append(rec)
    return out


FAMILY_MIN_WEBS = {"A": 0, "B": 0, "C": 1, "D": 2}
T_SPAR_MIN = 0.002  # STRUCTURAL_MIN_T_SPAR


def greedy_drop_to_usable(family: str, airfoil, survivors: list[float], layout) -> dict[str, Any]:
    """Drop the aft-most web until the cap allowance reaches the material minimum.

    Returns how many extra drops that costs beyond bare geometric feasibility,
    and whether the family's minimum web count blocks it first.
    """
    if not bool(layout.has_spar_caps) or not survivors:
        return {"applies": False}
    floor = FAMILY_MIN_WEBS.get(family, 0)
    xs = list(survivors)
    for extra in range(len(xs)):
        cap = cap_probe(family, airfoil, xs, layout)
        allowed = float(cap.get("cap_h_allowed_m", 0.0))
        if allowed >= T_SPAR_MIN:
            return {"applies": True, "extra_drops": extra, "n_webs": len(xs),
                    "cap_h_m": allowed, "blocked_by_floor": False}
        if len(xs) - 1 < max(floor, 1):
            return {"applies": True, "extra_drops": extra, "n_webs": len(xs),
                    "cap_h_m": allowed, "blocked_by_floor": True}
        xs = xs[:-1]  # drop the aft-most
    return {"applies": True, "extra_drops": len(survivors) - 1, "n_webs": 1,
            "cap_h_m": 0.0, "blocked_by_floor": True}


def cap_probe(family: str, airfoil, xs_resolved: list[float], layout) -> dict[str, Any]:
    """Cap band extent and clamped height, following build_spar_caps' family split."""
    if not bool(layout.has_spar_caps) or not xs_resolved:
        return {"applies": False}
    common = dict(y_search=None, phi_inner_ref=None)
    if family == "D":
        # Continuous band, first web to last, height keyed on the aft end.
        x_lo, x_hi = float(xs_resolved[0]), float(xs_resolved[-1])
        _, _, h_u, h_l, info = _cap_feasibility_clamp(
            airfoil, T_SKIN, x_lo, x_hi, T_SPAR, T_SPAR,
            aft_half_depth_rule=True, **common,
        )
        return {
            "applies": True, "mode": "continuous_aft_half_depth",
            "x_lo": x_lo, "x_hi": x_hi, "band_m": x_hi - x_lo,
            "cap_h_allowed_m": float(h_u), "requested_m": T_SPAR,
            "scaled": bool(info.get("scaled_heights")),
        }
    if family == "C":
        # One discrete cap per web; the common minimum governs them all.
        dhw = 0.05 * float(airfoil.chord)
        per = []
        for x in xs_resolved:
            _, _, h_u, _h_l, _i = _cap_feasibility_clamp(
                airfoil, T_SKIN, float(x) - dhw, float(x) + dhw, T_SPAR, T_SPAR,
                aft_half_depth_rule=False, **common,
            )
            per.append(float(h_u))
        return {
            "applies": True, "mode": "discrete_per_web", "n_caps": len(per),
            "cap_h_per_band_m": per,
            "cap_h_allowed_m": min(per) if per else 0.0, "requested_m": T_SPAR,
        }
    # Family B: band centred on max thickness, independent of the webs.
    from blade_geometry_model.section_geometry.layout.section_axes import (
        max_thickness_chord_x,
    )

    fhw = 0.15 * float(airfoil.chord)
    x_c = float(max_thickness_chord_x(airfoil))
    _, _, h_u, _h_l, _i = _cap_feasibility_clamp(
        airfoil, T_SKIN, x_c - fhw, x_c + fhw, T_SPAR, T_SPAR,
        aft_half_depth_rule=False, **common,
    )
    return {
        "applies": True, "mode": "web_independent_band",
        "x_lo": x_c - fhw, "x_hi": x_c + fhw, "band_m": 2 * fhw,
        "cap_h_allowed_m": float(h_u), "requested_m": T_SPAR,
    }


def report(layout_key: str, recs: list[dict[str, Any]]) -> None:
    print(f"\n{'=' * 100}\n{layout_key}\n{'=' * 100}")
    hdr = (f"{'stn':>4}{'z/L':>7}{'chord_mm':>10}{'nweb':>6}{'drop':>6}"
           f"{'maxwalk%c':>11}{'cap_h_mm':>10}"
           f"{'|':>3}{'RA_nweb':>9}{'RA_drop':>9}{'RA_cap_h':>10}")
    print(hdr)
    print("-" * len(hdr))
    for r in recs:
        if "error" in r:
            print(f"{r['idx']:>4}{r['zL']:>7.3f}   ERROR {r['error'][:60]}")
            continue
        cap = r.get("cap") or {}
        ra = r.get("root_anchored") or {}
        ra_cap = ra.get("cap") or {}
        print(
            f"{r['idx']:>4}{r['zL']:>7.3f}{r['chord'] * 1e3:>10.1f}"
            f"{r['n_webs_resolved']:>6}{r['n_dropped']:>6}"
            f"{r['max_walk_frac_chord'] * 100:>11.2f}"
            f"{cap.get('cap_h_allowed_m', float('nan')) * 1e3:>10.3f}"
            f"{'|':>3}{ra.get('n_survivors', -1):>9}{ra.get('n_dropped', -1):>9}"
            f"{ra_cap.get('cap_h_allowed_m', float('nan')) * 1e3:>10.3f}"
        )
    ok = [r for r in recs if "error" not in r]
    if not ok:
        return
    dropped = [r for r in ok if r["n_dropped"] > 0]
    walked = [r for r in ok if r["max_walk_frac_chord"] > 1e-9]
    short = [r for r in ok
             if (r.get("cap") or {}).get("cap_h_allowed_m", 1.0) < 0.99 * T_SPAR]
    print(f"\n  stations                          : {len(ok)}")
    print(f"  with >=1 web dropped              : {len(dropped)}")
    print(f"  with >=1 web walked LE-ward       : {len(walked)}")
    if walked:
        print(f"  worst LE-ward walk                : "
              f"{max(r['max_walk_frac_chord'] for r in walked) * 100:.2f}% of chord")
    print(f"  where cap height was clamped      : {len(short)} "
          f"(requested {T_SPAR * 1e3:.2f} mm)")
    if short:
        worst = min(short, key=lambda r: r["cap"]["cap_h_allowed_m"])
        print(f"  worst clamp                       : station {worst['idx']} "
              f"(z/L={worst['zL']:.3f}) -> "
              f"{worst['cap']['cap_h_allowed_m'] * 1e3:.4f} mm")

    ra_ok = [r for r in ok if r.get("root_anchored")]
    if ra_ok:
        ra_drop = [r for r in ra_ok if r["root_anchored"]["n_dropped"] > 0]
        print("  -- root-anchored counterfactual --")
        print(f"  stations losing >=1 web           : {len(ra_drop)}")
        if ra_drop:
            first = min(ra_drop, key=lambda r: r["idx"])
            print(f"  first drop at                     : station {first['idx']} "
                  f"(z/L={first['zL']:.3f}), "
                  f"{first['root_anchored']['n_survivors']} of "
                  f"{first['n_webs_nominal']} webs survive")
            counts = {r["idx"]: r["root_anchored"]["n_survivors"] for r in ra_ok}
            seq = [counts[i] for i in sorted(counts)]
            print(f"  survivor count root -> tip        : "
                  f"{seq[0]} -> {seq[-1]} (non-increasing: "
                  f"{all(b <= a for a, b in zip(seq, seq[1:]))})")
        ra_short = [r for r in ra_ok
                    if (r["root_anchored"].get("cap") or {}).get(
                        "cap_h_allowed_m", 1.0) < 0.99 * T_SPAR]
        print(f"  where cap height would clamp      : {len(ra_short)}")
        if ra_short:
            w = min(ra_short, key=lambda r: r["root_anchored"]["cap"]["cap_h_allowed_m"])
            print(f"  worst clamp                       : station {w['idx']} -> "
                  f"{w['root_anchored']['cap']['cap_h_allowed_m'] * 1e3:.4f} mm")

    g_ok = [r for r in ok if (r.get("greedy") or {}).get("applies")]
    if g_ok:
        blocked = [r for r in g_ok if r["greedy"]["blocked_by_floor"]]
        extra = [r for r in g_ok if r["greedy"]["extra_drops"] > 0]
        print(f"  -- greedy drop to cap >= {T_SPAR_MIN * 1e3:.0f} mm --")
        print(f"  stations needing extra drops      : {len(extra)}")
        print(f"  stations blocked by family floor  : {len(blocked)}")
        if extra:
            print(f"  max extra drops at one station    : "
                  f"{max(r['greedy']['extra_drops'] for r in extra)}")
        seq = [r["greedy"]["n_webs"] for r in sorted(g_ok, key=lambda r: r["idx"])]
        print(f"  final web count root -> tip       : {seq[0]} -> {seq[-1]} "
              f"(non-increasing: {all(b <= a for a, b in zip(seq, seq[1:]))})")
        if blocked:
            w = min(blocked, key=lambda r: r["greedy"]["cap_h_m"])
            print(f"  worst blocked station             : {w['idx']} "
                  f"({w['greedy']['n_webs']} webs, cap "
                  f"{w['greedy']['cap_h_m'] * 1e3:.4f} mm)")


def main() -> int:
    stations = build_stations()
    print(f"reference rotor {ROTOR}: {len(stations)} stations, "
          f"chord {stations[0]['chord'] * 1e3:.1f} -> {stations[-1]['chord'] * 1e3:.1f} mm")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    everything = {}
    for key in LAYOUTS:
        recs = probe_layout(key, stations)
        report(key, recs)
        everything[key] = recs
    path = OUT_DIR / "spanwise_web_anchors.json"
    path.write_text(json.dumps(everything, indent=1, default=float), encoding="utf-8")
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
