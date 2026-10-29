# Frames inventory — coordinate-frame safety baseline

**Programme:** Coordinate-Frame Safety (WP0 baseline)  
**Date:** 2026-09-23  
**Authority:** this file is the programme status source. Semantics and transforms
remain in [`README.md`](README.md) / [`section_frames/`](section_frames/) /
[`shell_frames/`](shell_frames/). Cross-cutting MUST rules:
[`docs/conventions/api/API_STANDARDS.md`](../../../../../docs/conventions/api/API_STANDARDS.md)
and `ARCH-MUST-016` in
[`MODULE_ARCHITECTURE_STANDARDS.md`](../../../../../docs/conventions/api/MODULE_ARCHITECTURE_STANDARDS.md).

**Tagging policy:** validated `FrameTag` is required at public / cross-package
boundaries carrying frame-dependent geometry, resultants, stresses, stiffnesses,
or tensors. Documented fixed-frame internal arrays may remain untagged.

**Frozen exceptions (not FrameTag members):**

| Exception | Rationale |
|-----------|-----------|
| NACA / SDF / chord authoring space | Pre-handoff geometry; label as authoring space, not `"geometry"` FrameTag |
| C / isoparametric ξ,η | Element kinematics |
| FE assembly index space | DOF numbering, not structural B/S |
| External-tool adapters | Cotide / sectionproperties bridges |
| Display-only plot axis invert | Does not change stored `(X_*, Y_*)` |
| WP4.5 aliases (`y`/`z`, `FRAME_SECTION_YZ`) | Dated compatibility; see [`WP45_SECTION_FRAME_RENAME_MAP.md`](../../../../../docs/plans/long_term_structure/reports/WP45_SECTION_FRAME_RENAME_MAP.md) |

Classification: `frame-dependent` | `fixed-frame internal` | `not frame-dependent`.

Status: `compliant` | `remediation-required` | `N/A` | `dated exception`.

---

## Package classification

### `blade_structure_precompute/module`

| Package | Class | Status | Notes / WP |
|---------|-------|--------|------------|
| `section_selection` | frame-dependent | compliant | FRAME-002/006/007 closed 2026-09-23 |
| `stiffness_selection` | frame-dependent | compliant | Uses `CoordinateFrameTransforms`; B via naming |
| `feasibility` | frame-dependent | compliant | FRAME-001 closed 2026-09-23 |
| `internal_force` | fixed-frame internal | compliant | B via `_XB`/`_YB`/`_ZB` suffixes; no local transforms |
| `utilities/stress` | frame-dependent | compliant | FRAME-003 closed — B touch lever arms via `rotate_S_to_B` |
| `utilities/effective_width_analyzer` | frame-dependent | compliant | Delegates B→S moments |
| `utilities/required_section_depth` | frame-dependent | compliant | Documents B moment; uses decompose helper |
| `utilities/effective_width_lookup` | not frame-dependent | N/A | Tabular lookup |
| `utilities/spanwise_schedule` | not frame-dependent | N/A | Span schedule only |
| `utilities/core_budget` | not frame-dependent | N/A | Scalar budget |

### `blade_structure_utilities`

| Package | Class | Status | Notes / WP |
|---------|-------|--------|------------|
| `geometry_utilities/coordinate_frames` | frame-dependent | compliant | SoT |
| `geometry_utilities/section_contract` | frame-dependent | compliant | `BeamResultantVector.frame`; B midlines |
| `geometry_utilities/section_geometry` | frame-dependent | compliant | Legacy `rotate_chord_to_blade` removed; midlines use `chord_to_blade` |
| `geometry_utilities/section_properties` | frame-dependent | compliant | Authoring-space vs FrameTag via `frame_labels.normalize_section_properties_frame` |
| `geometry_utilities/section_models` | fixed-frame internal | compliant | Authoring models |
| `geometry_utilities/scale` | frame-dependent | compliant | Delegates transforms |
| `geometry_utilities/root_blend_model` | frame-dependent | dated exception | `frame: str` literals; review 2026-12-31 |
| `mesh_utilities/section_mesh` | frame-dependent | compliant | Canon B handoff + tests |
| `mesh_utilities/area_2d` | not frame-dependent | N/A | Kernel |
| `mesh_utilities/line_1d` | not frame-dependent | N/A | Kernel |
| `load_utilities/section_loads` | frame-dependent | compliant | Tagged resultants |
| `load_utilities/stress` | frame-dependent | compliant | Contour EB; S convention documented |
| `load_utilities/internal_force` | fixed-frame internal | compliant | B naming |
| `fem_utilities/section_shell_fea` | frame-dependent | compliant | Strongest canon tests |
| `material_utilities/laminate` | frame-dependent | compliant | E↔M via canon |
| `material_utilities/orthotropic` | frame-dependent | compliant | Legacy `rotation_matrix_theta` removed; use `ply_stress_transform_E_to_M` |
| `material_utilities/section_bind` | fixed-frame internal | compliant | Strip materials |
| `material_utilities/section_gauge` | fixed-frame internal | compliant | Count-native gauges |
| `material_utilities/section_layup` | fixed-frame internal | dated exception | Local shear rotate; review WP3/WP5 |
| `material_utilities/section_schedule` | not frame-dependent | N/A | Span reconcile |
| `material_utilities/section_sizing` | fixed-frame internal | compliant | S-authoring depth |
| `material_utilities/failure` | fixed-frame internal | compliant | FI in M via callers |
| `material_utilities/isotropic` | not frame-dependent | N/A | |
| `material_utilities/layup` | fixed-frame internal | compliant | Stack types |
| `integration_utilities/*` | fixed-frame internal | compliant | S by API name |
| `fatigue_utilities/*` | not frame-dependent | N/A | Stress histories tagged upstream |
| `parallel_utilities` | not frame-dependent | N/A | |
| `logging_utilities` | not frame-dependent | N/A | |

---

## Finding register

| ID | Severity | Owner | Evidence | Target WP | Status |
|----|----------|-------|----------|-----------|--------|
| FRAME-001 | critical | feasibility | `tc_iteration.py` inline S→B omits `Ixy` | WP2 | **closed** |
| FRAME-002 | high | section_selection | silent S-as-B in evaluation / property packing | WP2 | **closed** |
| FRAME-003 | high | utilities/stress | `compute_stress_at_touch_points_B` uses S touch X as B | WP2 | **closed** |
| FRAME-004 | high | section_geometry | exported unused `rotate_chord_to_blade` | WP3 | **closed** (removed 2026-09-23) |
| FRAME-005 | medium | orthotropic | public `rotation_matrix_theta` duplicates/misnames Tσ | WP3 | **closed** (removed 2026-09-23) |
| FRAME-006 | medium | section_selection | `SectionAppliedLoadData` / row loads untagged | WP4 | **closed** |
| FRAME-007 | medium | section_selection | `SectionFramePropertyData` typed but unvalidated | WP4 | **closed** |
| FRAME-008 | medium | section_properties | `frame="geometry"` ad-hoc string | WP4 | **closed** (authoring-space vocabulary) |
| FRAME-009 | low | conventions | no ARCH-MUST / template DoD for frames | WP1 | **closed** |
| FRAME-010 | low | coordinate_frames | missing `FRAMES_INVENTORY.md` (this file) | WP0 | **closed** |
| FRAME-011 | low | enforcement | no repo-wide transform-delegation AST gate | WP6 | **closed** |

---

## Baseline tests (WP0)

Existing focused coverage before remediation:

| Suite | Path |
|-------|------|
| Frame tags | `geometry_utilities/coordinate_frames/tests/test_frame_tags.py` |
| Section plane | `mesh_utilities/section_mesh/tests/test_section_plane_canon.py` |
| FE canon | `fem_utilities/section_shell_fea/tests/test_section_frame_canon.py` |
| Beam resultant tags | `geometry_utilities/section_contract/tests/test_beam_resultant_frames.py` |
| Ply E↔M | `geometry_utilities/coordinate_frames/tests/test_ply_frame_transforms.py` |
| Campaign S/B packing | `blade_structure_precompute/.../test_campaign_section_property_frames.py` |

Run (non-slow):

```bash
pytest -q -m "not slow" src/blade_structure_utilities/geometry_utilities/coordinate_frames/tests src/blade_structure_utilities/mesh_utilities/section_mesh/tests/test_section_plane_canon.py src/blade_structure_utilities/fem_utilities/section_shell_fea/tests/test_section_frame_canon.py src/blade_structure_utilities/geometry_utilities/section_contract/tests/test_beam_resultant_frames.py
```

---

## Status updates

Update this table as each WP closes a finding. Do not delete rows.
