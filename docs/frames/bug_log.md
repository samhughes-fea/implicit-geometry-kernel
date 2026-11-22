# bug_log — coordinate_frames

**Audience: internal.** Maintainer notes; not part of the library manual.

↑ [`coordinate_frames`](README.md)

## Role

Package debug log. Identifiers are `CF-nnn`. Cross-link
[`section/bug_log.md`](../section/bug_log.md) (`SG-`) and
[`blade/bug_log.md`](../blade/bug_log.md) (`BG-`) when a defect spans packages.

## Index

| ID | Title | Status |
|---|---|---|
| [CF-001](#cf-001) | No shared planar rotation harness for B-consuming packages | **Fixed** (2026-09-25) |
| [CF-002](#cf-002) | `transform_stiffness_B_to_Sprime` invented 1e-6 / 1e-12 EI | **Fixed** (2026-09-26) |
| [CF-003](#cf-003) | `transform_second_moments_S_to_B` wrong Ixy sign on `Ixx_B` | **Fixed** (2026-11-05) |
| [CF-004](#cf-004) | `transform_moments_S_to_B` / `B_to_S` rotate by the opposite twist to points and forces | **Fixed** (2026-11-22) |

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

---

<a id="cf-004"></a>
## CF-004 — Moment vectors rotate by the opposite twist to points and forces

**Status:** Fixed — 2026-11-22 (found the same day while deriving [R1](../r/R1_frames.md))
**Area:** `frames.tensors.moments` (`transform_moments_B_to_S`,
`transform_moments_S_to_B`)

### Symptom
For a non-zero twist, `transform_moments_S_to_B` returns a different moment from the
one obtained by moving the loaded points to B with `rotate_S_to_B` and taking the
moments there. Points (`rotate_S_to_B`), forces (`transform_forces_S_to_B`) and
second moments (`transform_second_moments_S_to_B`) all agree with each other; the
moment transform agrees with them only if the sign of the twist is reversed.

### Evidence
Fifty axial forces \(F_z\) at random points in S, twist 20°. The in-plane moments
\(M_X=\sum Y F_z\), \(M_Y=-\sum X F_z\) computed in B after `rotate_S_to_B`:
\((M_{X_B}, M_{Y_B}) = (-1.2027, -7.6630)\). `transform_moments_S_to_B` with
`twist_deg=20`: \((4.0044, -6.6432)\). With `twist_deg=-20`: \((-1.2027, -7.6630)\),
equal to 1e-15. The derivation in R1, (R1.15)–(R1.16), gives
\(\mathbf{M}_B=\mathbf{Q}(\theta)\mathbf{M}_S\), the same rotation as points.

### Root cause (diagnosis)
The moment functions take their components in the order \((M_Y, M_X)\) but apply the
matrix written for \((X, Y)\) order. Swapping the components reverses the rotation,
\(\mathbf{P}\mathbf{Q}(\theta)\mathbf{P}=\mathbf{Q}(-\theta)\) (R1.10), the same trap
documented on `rotate_S_to_B_swapped_order`. The unit test
`tests/frames/unit/test_coordinate_force_transforms.py::test_matches_moment_transform`
asserts the force and moment transforms agree after passing `M_YB=V_XB`, which pins
the swapped behaviour rather than checking it against the point rotation.

### Fix
`transform_moments_B_to_S` now forms \((M_{X_S}, M_{Y_S}) = \mathbf{Q}(\theta)^{\mathsf{T}}(M_{X_B}, M_{Y_B})\)
with the components written out by name, so the sine terms change sign;
`transform_moments_S_to_B` inverts it. `transform_forces_B_to_S` no longer relies on
the old argument order: it passes \((V_Y, V_X)\) to the moment transform. Signatures
and argument order are unchanged; results change for any non-zero twist (breaking,
recorded in the changelog). No public sibling package calls these functions (checked
2026-11-22). Decided by the maintainer on 2026-11-22.

### Verification
`tests/frames/unit/test_moments_rotate_with_points.py`: moments taken from points
rotated by `rotate_S_to_B` equal `transform_moments_S_to_B` to 1e-12 at eight twists
from −90° to 90°, and B → S inverts it; a unit flapwise moment, force and point
rotate to the same vector. The force–moment agreement test now passes components in
the same \((X, Y)\) order.

