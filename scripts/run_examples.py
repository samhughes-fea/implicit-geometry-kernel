"""Run every example and demo script headless and report the ones that fail.

    python scripts/run_examples.py            # all scripts (blade/09 --smoke)
    python scripts/run_examples.py --full     # blade/09 every catalogue key
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
# Extra argv for the CI/default run. --full drops these so maintainers can
# regenerate the complete blade/09 gallery.
SMOKE_ARGS = {
    "blade/09_blade_concepts/blade_concepts.py": ["--smoke"],
}
# Wall-clock seconds when --full runs the complete 09 catalogue.
SLOW = {
    "blade/09_blade_concepts/blade_concepts.py": 720.0,
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


def run(path: Path, timeout: float, extra: list[str] | None = None, full: bool = False) -> tuple[Path, int, str]:
    env = dict(os.environ, MPLBACKEND="Agg")
    rel = path.relative_to(EXAMPLES).as_posix()
    if full:
        timeout = max(float(timeout), float(SLOW.get(rel, 0.0)))
    extra = list(extra or ())
    try:
        proc = subprocess.run([sys.executable, str(path), *extra], cwd=path.parent, env=env, capture_output=True,
                              text=True, timeout=timeout)
        if proc.returncode:
            err = (proc.stderr or proc.stdout or "").strip()
            return path, proc.returncode, err[-4000:] if err else "no output"
        return path, 0, ""
    except subprocess.TimeoutExpired:
        return path, -1, f"timed out after {timeout:g} s"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--timeout", type=float, default=480.0, help="seconds per script")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    parser.add_argument("--full", action="store_true",
                        help="run blade/09 for every catalogue key (slow)")
    args = parser.parse_args()
    paths = scripts()

    def _job(path: Path) -> tuple[Path, int, str]:
        rel = path.relative_to(EXAMPLES).as_posix()
        extra = [] if args.full else SMOKE_ARGS.get(rel, [])
        return run(path, args.timeout, extra, full=args.full)

    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(_job, paths))
    failed = [(p, code, msg) for p, code, msg in results if code != 0]
    for p, code, msg in failed:
        print(f"FAIL {p.relative_to(ROOT)} (exit {code})")
        print(msg)
    print(f"{len(paths) - len(failed)} of {len(paths)} scripts ran ({len(SKIP)} skipped)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
