"""Test configuration helpers for environments without pytest-cov."""

from __future__ import annotations

import sys
from pathlib import Path


# Ensure the repository root is importable without ``pip install -e .``
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:  # pragma: no cover - defensive path injection
    sys.path.insert(0, str(ROOT))


try:  # pragma: no cover - best-effort detection
    import pytest_cov  # type: ignore
except ImportError:  # pragma: no cover - executed when plugin missing
    def pytest_addoption(parser):
        parser.addoption("--cov", action="append", default=[], help="No-op without pytest-cov")
        parser.addoption("--cov-report", action="append", default=[], help="No-op without pytest-cov")
