"""Diff two example-15 figure folders (health_counts.csv + failures.csv).

Classify-only: prints count deltas and element ids that appeared or left.
No pass/fail gate.

    python examples/blade/15_shell_mesh/compare_studies.py \\
        figures/nrel_5mw figures/nrel_5mw_after
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from implicit_geometry_kernel.mesh import diff_study_folders  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path, help="Earlier study folder")
    parser.add_argument("after", type=Path, help="Later study folder")
    args = parser.parse_args(argv)
    report = diff_study_folders(args.before, args.after)
    print(f"compare {args.before} → {args.after}")
    if not report["counts"]:
        print("  (no health_counts rows in either folder)")
    for row in report["counts"]:
        sign = "+" if row["delta"] > 0 else ""
        print(
            f"  {row['code']}: {row['before']} → {row['after']} "
            f"({sign}{row['delta']})"
        )
    new_ids = report["elements"]["new"]
    gone_ids = report["elements"]["gone"]
    print(f"  elements new: {new_ids[:20]}{'…' if len(new_ids) > 20 else ''}")
    print(f"  elements gone: {gone_ids[:20]}{'…' if len(gone_ids) > 20 else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
