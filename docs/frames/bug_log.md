# bug_log — coordinate_frames

**Audience: internal.** Maintainer notes; not part of the library manual.

↑ [`coordinate_frames`](README.md)

## Role

Package debug log. Identifiers are `CF-nnn`. Cross-link repo-root `BL-nnn` and
sibling package logs (`SG-`, `SL-`, …) when a defect spans packages.

## Index

| ID | Title | Status |
|---|---|---|
| [CF-001](#cf-001) | No shared planar rotation harness for B-consuming packages | **Fixed** (2026-09-25) |
| [CF-002](#cf-002) | `transform_stiffness_B_to_Sprime` invented 1e-6 / 1e-12 EI | **Fixed** (2026-09-26) |
| [CF-003](#cf-003) | `transform_second_moments_S_to_B` wrong Ixy sign on `Ixx_B` | **Fixed** (2026-11-05) |

---

<a id="cf-001"></a>
## CF-001 — No shared planar rotation harness for B-consuming packages

**Status:** Fixed — 2026-09-25
**Area:** `coordinate_frames/testing/invariance.py`

### Symptom
B-consuming packages each invented (or skipped) rotation checks. Three Stage 0
defects (SG-014, SG-015, SL-008) were chord-first code on B-frame input or an
adapter that mirrored without converting handedness.

### Evidence
Exactness programme item 2 survey: few package-local tests; no shared property
harness that rotates a section and requires topology / stiffness / loads to
transform correctly.

### Root cause
`coordinate_frames` owned the transforms but not a reusable invariance test
surface for consumers.

### Fix
`testing/invariance.py` exports `rotate_midlines_in_b`,
`assert_planar_pierce_rotates` and `assert_scalar_invariant_under_rotation`, all
gated on `FRAME_B` / `require_frame_tag`. Pinned by
`tests/test_frame_invariance_harness.py` on the thin-wall shear centre and `GJ`.

### Verification
`pytest …/coordinate_frames/tests/test_frame_invariance_harness.py` green.

---

<a id="cf-002"></a>
## CF-002 — Silent 1e-6 / 1e-12 EI invent in B→S′ transform

**Status:** Fixed — 2026-09-26
**Area:** `section_frames/stiffness.py`

### Symptom
`except Exception` on shoelace / polygon failure wrote `EI_edge_B = 1e-6` or
`I = 1e-12/E`; `n_points < 3` skipped or used pivot EI; negative centroid I was
clamped to `1e-12`. Missing `EI_edge_B` used `Ixx_S` as if it were `Ixx_B`.

### Fix
Raise with CF-002. Derive missing edge EI via `transform_second_moments_S_to_B`.
No 1e-12 clamp.

### Verification
`tests/test_cf002_sp001_exactness.py`

---

<a id="cf-003"></a>
## CF-003 — Wrong Ixy sign in `Ixx_B` / `Iyy_B` under twist

**Status:** Fixed — 2026-11-05
**Area:** `frames.tensors.inertia` (`transform_second_moments_S_to_B` /
`transform_second_moments_B_to_S`)

### Symptom
Rotating second moments S→B did not match integrating the outline after
`rotate_S_to_B`. `Ixx_B + Iyy_B` drifted with twist (trace not invariant).

### Evidence
`tests/frames/unit/test_second_moments_rotate_with_points.py` (added with the
fix): five twists on a NACA 2418 outline; tensor transform must agree with
shoelace on the rotated polygon, and `Ixx + Iyy` must match the S-frame sum.
Logged in `CHANGELOG.md` 0.4.0; missing from this package log until now.

### Root cause
Point rotation is \(X_B = c X_S + s Y_S\), \(Y_B = -s X_S + c Y_S\). Substituting
into \(I_{xx}=\int Y^2\) gives a **minus** \(2 I_{xy} s c\) term. The transform
used **plus** \(2 I_{xy} \sin\cos\) on `Ixx_B` (and the opposite error on
`Iyy_B`), so the product-moment contribution had the wrong sign.

### Fix
`Ixx_B = Ixx_S cos² + Iyy_S sin² − 2 Ixy_S sin cos` (and the matching `Iyy_B` /
`Ixy_B` forms) in `inertia.py`. Inverse `B→S` uses the same convention.
Commit `52aa320`.

### Verification
`pytest tests/frames/unit/test_second_moments_rotate_with_points.py` green;
`Ixx_B + Iyy_B ≈ Ixx_S + Iyy_S` at every twist in that suite.

