"""Run a downstream package's tests against this checkout of the kernel.

    python scripts/check_consumer.py PATH [PYTEST ARGS ...]

PATH is a local checkout of any package that depends on the kernel. The script
installs that package without its dependencies (so the kernel installed here, this
checkout, is the one it imports) and runs its tests. Extra arguments go to pytest,
for example ``-k geometry`` or a list of test files. Exit status is pytest's.

Use it before merging a change that could affect callers. The kernel names no
consumer; which packages to check is up to whoever runs it.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main(argv: list[str]) -> int:
    if not argv or argv[0] in {"-h", "--help"}:
        print(__doc__)
        return 0 if argv else 2
    consumer = Path(argv[0]).resolve()
    if not (consumer / "pyproject.toml").is_file() and not (consumer / "setup.py").is_file():
        print(f"{consumer} is not a Python package checkout")
        return 2
    py = sys.executable
    subprocess.run([py, "-m", "pip", "install", "-q", "-e", str(ROOT)], check=True)
    subprocess.run([py, "-m", "pip", "install", "-q", "--no-deps", "-e", str(consumer)], check=True)
    kernel = subprocess.run(
        [py, "-c", "import implicit_geometry_kernel as k; print(k.__file__)"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    if not Path(kernel).resolve().is_relative_to(ROOT):
        print(f"the consumer resolved the kernel from {kernel}, not this checkout")
        return 2
    print(f"kernel: {ROOT}\nconsumer: {consumer}")
    return subprocess.run([py, "-m", "pytest", "-q", "-p", "no:cacheprovider", *argv[1:]], cwd=consumer).returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
