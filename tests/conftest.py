"""Shared pytest hooks."""

from __future__ import annotations

import pytest


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    # Tests that need the private reference-blades catalogue carry a skipif whose
    # reason names it. Mark them ``reference`` so ``pytest -m reference`` selects them.
    for item in items:
        for mark in item.iter_markers("skipif"):
            if "catalogue" in str(mark.kwargs.get("reason", "")):
                item.add_marker(pytest.mark.reference)
                break
