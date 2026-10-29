# bug_log — coordinate_frames

↑ [`coordinate_frames`](../README.md)

## Role

Package debug log. Identifiers are `CF-nnn`. Cross-link repo-root `BL-nnn` and
sibling package logs (`SG-`, `SL-`, …) when a defect spans packages.

## Index

| ID | Title | Status |
|---|---|---|
| [CF-001](#cf-001) | No shared planar rotation harness for B-consuming packages | **Fixed** (2026-09-25) |
| [CF-002](#cf-002) | `transform_stiffness_B_to_Sprime` invented 1e-6 / 1e-12 EI | **Fixed** (2026-09-26) |

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

