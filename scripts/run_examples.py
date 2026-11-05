"""Run every example and demo script headless and report the ones that fail.

    python scripts/run_examples.py            # all scripts
    python scripts/run_examples.py --timeout 300

Scripts listed in SKIP are known to be broken or too slow for CI; each says why.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"

SKIP = {
    # the hybrid props report no longer exposes a residual quadtree mesh
    "demos/section/le_te_insert_residual_demo.py",
    # renders every layout at full resolution; takes several minutes
    "demos/section/quadtree_union_mesh.py",
}
HELPERS = {"_ensure_repo_path.py", "_plot_debug.py", "_taxonomy_cases.py", "cell_wall_loops_case_payload.py",
           "feature_flythrough_specs.py"}


def scripts() -> list[Path]:
    out = []
    for path in sorted(EXAMPLES.rglob("*.py")):
        rel = path.relative_to(EXAMPLES).as_posix()
        if path.name in HELPERS or rel in SKIP or "output" in path.parts or "_shared" in path.parts:
            continue
        out.append(path)
    return out


def run(path: Path, timeout: float) -> tuple[Path, int, str]:
    env = dict(os.environ, MPLBACKEND="Agg")
    try:
        proc = subprocess.run([sys.executable, str(path)], cwd=path.parent, env=env, capture_output=True,
                              text=True, timeout=timeout)
        tail = (proc.stderr or proc.stdout).strip().splitlines()[-1:] if proc.returncode else []
        return path, proc.returncode, tail[0] if tail else ""
    except subprocess.TimeoutExpired:
        return path, -1, f"timed out after {timeout:g} s"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--timeout", type=float, default=240.0, help="seconds per script")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    args = parser.parse_args()
    paths = scripts()
    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(lambda p: run(p, args.timeout), paths))
    failed = [(p, code, msg) for p, code, msg in results if code != 0]
    for p, code, msg in failed:
        print(f"FAIL {p.relative_to(ROOT)}: {msg}")
    print(f"{len(paths) - len(failed)} of {len(paths)} scripts ran ({len(SKIP)} skipped)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
