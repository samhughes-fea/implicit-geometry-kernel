"""Scenario and profile definitions for section-geometry fly-throughs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

FrameSource = Literal["spanwise", "synthetic"]
PanelLayout = Literal["single", "dual", "tri"]


@dataclass(frozen=True)
class RenderProfile:
    name: str
    width_px: int
    frame_delay_cs: int
    dpi: int
    hold_head: int
    hold_tail: int
    grid_nx: int
    grid_ny: int


@dataclass(frozen=True)
class ScenarioSpec:
    name: str
    output_name: str
    title: str
    frame_source: FrameSource
    panel_layout: PanelLayout
    description: str
    frame_count: int
    sweep_params: dict[str, Any]
    metric_overlays: tuple[str, ...]


PROFILES: dict[str, RenderProfile] = {
    "preview": RenderProfile(
        name="preview",
        width_px=1024,
        frame_delay_cs=12,
        dpi=130,
        hold_head=3,
        hold_tail=4,
        grid_nx=280,
        grid_ny=120,
    ),
    "balanced": RenderProfile(
        name="balanced",
        width_px=1280,
        frame_delay_cs=9,
        dpi=155,
        hold_head=5,
        hold_tail=7,
        grid_nx=420,
        grid_ny=170,
    ),
    "high": RenderProfile(
        name="high",
        width_px=1600,
        frame_delay_cs=7,
        dpi=200,
        hold_head=6,
        hold_tail=8,
        grid_nx=560,
        grid_ny=240,
    ),
}

DEFAULT_PROFILE = "balanced"
DEFAULT_ROTOR = "new_3.2m_rotor"

SCENARIOS: dict[str, ScenarioSpec] = {
    "anchoring_cascade": ScenarioSpec(
        name="anchoring_cascade",
        output_name="feature_01_spanwise_anchoring_cascade.gif",
        title="Spanwise Anchoring Cascade",
        frame_source="spanwise",
        panel_layout="dual",
        description="Station-relative vs root-anchored web behavior along span.",
        frame_count=64,
        sweep_params={"driver": "station_index", "domain": "new_3.2m_rotor span"},
        metric_overlays=("station", "zL", "chord", "survivors", "drops"),
    ),
    "key_cascade_spanwise": ScenarioSpec(
        name="key_cascade_spanwise",
        output_name="feature_02_structural_key_cascade_spanwise.gif",
        title="Structural Key Cascade (Spanwise)",
        frame_source="spanwise",
        panel_layout="tri",
        description="Spanwise comparison of 5B-F, 5C-CN, and 5D-F section systems.",
        frame_count=64,
        sweep_params={"driver": "station_index", "keys": ["5B-F", "5C-CN", "5D-F"]},
        metric_overlays=("station", "zL", "root_survivors", "cap_allowance_mm"),
    ),
    "topology_twist": ScenarioSpec(
        name="topology_twist",
        output_name="feature_03_topology_under_twist.gif",
        title="Topology Under Twist",
        frame_source="synthetic",
        panel_layout="dual",
        description="Raw vs snapped strip topology through a twist sweep.",
        frame_count=31,
        sweep_params={"twist_start_deg": 0.0, "twist_end_deg": 30.0},
        metric_overlays=("twist_deg", "raw_gap_m", "snapped_gap_m"),
    ),
    "orientation_comparator": ScenarioSpec(
        name="orientation_comparator",
        output_name="feature_04_orientation_comparator.gif",
        title="Orientation Comparator (CN vs F)",
        frame_source="spanwise",
        panel_layout="dual",
        description="Matched station fly-through for chord-normal vs flapwise webs.",
        frame_count=64,
        sweep_params={"driver": "station_index", "pair": ["3D-CN", "3D-F"]},
        metric_overlays=("station", "zL", "orientation_pair"),
    ),
    "insert_dynamics": ScenarioSpec(
        name="insert_dynamics",
        output_name="feature_05_insert_dynamics_dual_rep.gif",
        title="Insert Dynamics (Dual-Rep)",
        frame_source="synthetic",
        panel_layout="dual",
        description="LE/TE insert sweep with hybrid residual-path visualization.",
        frame_count=18,
        sweep_params={
            "leading_edge_insert": [0.07, 0.15],
            "trailing_edge_insert": [0.72, 0.84],
        },
        metric_overlays=("le_insert_m", "te_insert_m", "residual_labels"),
    ),
    "quadtree_refinement": ScenarioSpec(
        name="quadtree_refinement",
        output_name="feature_06_quadtree_refinement_flythrough.gif",
        title="Quadtree Refinement Fly-Through",
        frame_source="synthetic",
        panel_layout="dual",
        description="Depth 4→9→4 refinement sweep with geometry and mesh panels.",
        frame_count=11,
        sweep_params={"depth_sequence": [4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4]},
        metric_overlays=("depth", "leaf_count"),
    ),
}

SCENARIO_ORDER = tuple(SCENARIOS.keys())


__all__ = [
    "DEFAULT_PROFILE",
    "DEFAULT_ROTOR",
    "PROFILES",
    "SCENARIOS",
    "SCENARIO_ORDER",
    "FrameSource",
    "PanelLayout",
    "RenderProfile",
    "ScenarioSpec",
]
