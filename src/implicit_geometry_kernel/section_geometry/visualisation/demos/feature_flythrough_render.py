"""Shared rendering helpers for section-geometry fly-through GIFs."""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from feature_flythrough_specs import RenderProfile


@dataclass(frozen=True)
class ScenarioPaths:
    scenario_name: str
    frame_dir: Path
    gif_path: Path
    metadata_path: Path


def prepare_scenario_paths(out_dir: Path, scenario_name: str, output_name: str) -> ScenarioPaths:
    """Create deterministic paths for one scenario output."""
    out_dir.mkdir(parents=True, exist_ok=True)
    frames_root = out_dir / "frames"
    frame_dir = frames_root / scenario_name
    frame_dir.mkdir(parents=True, exist_ok=True)
    return ScenarioPaths(
        scenario_name=scenario_name,
        frame_dir=frame_dir,
        gif_path=out_dir / output_name,
        metadata_path=out_dir / f"{scenario_name}.json",
    )


def reset_frame_dir(path: Path) -> None:
    """Delete stale frame files while preserving folder location."""
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def select_frame_indices(total: int, max_frames: int | None) -> list[int]:
    """Select deterministic frame indices, preserving endpoints."""
    if total <= 0:
        return []
    if max_frames is None or max_frames >= total:
        return list(range(total))
    if max_frames <= 1:
        return [0]

    raw = np.linspace(0, total - 1, max_frames)
    selected: list[int] = []
    for v in raw:
        idx = int(round(float(v)))
        if idx not in selected:
            selected.append(idx)

    if selected[0] != 0:
        selected.insert(0, 0)
    if selected[-1] != total - 1:
        selected.append(total - 1)

    if len(selected) > max_frames:
        selected = selected[:max_frames]
    while len(selected) < max_frames:
        for idx in range(total):
            if idx not in selected:
                selected.append(idx)
                if len(selected) == max_frames:
                    break
    return sorted(selected)


def frame_path(frame_dir: Path, frame_idx: int) -> Path:
    return frame_dir / f"frame_{frame_idx:04d}.png"


def apply_holds(
    frames: list[Path],
    *,
    head: int,
    tail: int,
) -> list[Path]:
    if not frames:
        return []
    return [frames[0]] * int(head) + frames + [frames[-1]] * int(tail)


def encode_gif(
    frame_paths: list[Path],
    out_path: Path,
    profile: RenderProfile,
) -> None:
    """Encode frames into an optimized looping GIF via ImageMagick."""
    if not frame_paths:
        raise ValueError("No frame paths supplied")

    list_file = out_path.with_suffix(".frames.txt")
    list_file.write_text("\n".join(str(p) for p in frame_paths), encoding="utf-8")
    try:
        cmd = [
            "magick",
            f"@{list_file}",
            "-coalesce",
            "-resize",
            f"{profile.width_px}x",
            "-colors",
            "256",
            "-layers",
            "Optimize",
            "-loop",
            "0",
            "-delay",
            str(profile.frame_delay_cs),
            str(out_path),
        ]
        subprocess.run(cmd, check=True)
    finally:
        list_file.unlink(missing_ok=True)


def validate_gif(path: Path) -> dict[str, Any]:
    """Return basic validation metrics for a GIF."""
    if not path.exists() or path.stat().st_size <= 0:
        raise RuntimeError(f"GIF missing or empty: {path}")

    identify = subprocess.run(
        ["magick", "identify", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    frame_count = sum(1 for line in identify.stdout.splitlines() if line.strip())

    first = subprocess.run(
        ["magick", "identify", "-format", "%w %h", f"{path}[0]"],
        check=True,
        capture_output=True,
        text=True,
    )
    width_s, height_s = first.stdout.strip().split()

    verbose = subprocess.run(
        ["magick", "identify", "-verbose", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    has_loop = "Iterations: 0" in verbose.stdout

    return {
        "path": str(path),
        "frame_count": int(frame_count),
        "width_px": int(width_s),
        "height_px": int(height_s),
        "infinite_loop": bool(has_loop),
        "size_bytes": int(path.stat().st_size),
    }


def write_metadata(path: Path, payload: dict[str, Any]) -> None:
    """Write deterministic scenario metadata JSON."""
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


__all__ = [
    "ScenarioPaths",
    "apply_holds",
    "encode_gif",
    "frame_path",
    "prepare_scenario_paths",
    "reset_frame_dir",
    "select_frame_indices",
    "validate_gif",
    "write_metadata",
]
