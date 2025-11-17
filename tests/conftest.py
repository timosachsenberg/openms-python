"""Test configuration helpers for environments without pytest-cov."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

try:  # pragma: no cover - best-effort detection
    import pytest_cov  # type: ignore
except ImportError:  # pragma: no cover - executed when plugin missing
    def pytest_addoption(parser):
        parser.addoption("--cov", action="append", default=[], help="No-op without pytest-cov")
        parser.addoption("--cov-report", action="append", default=[], help="No-op without pytest-cov")
