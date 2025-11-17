"""Shared helpers for OpenMS file IO handling."""
from __future__ import annotations

from pathlib import Path
from typing import FrozenSet, Iterable, Union


PathLike = Union[str, Path]

MS_EXPERIMENT_EXTENSIONS: FrozenSet[str] = frozenset({".mzml", ".mzxml", ".mzdata"})
FEATURE_MAP_EXTENSIONS: FrozenSet[str] = frozenset({".featurexml"})
CONSENSUS_MAP_EXTENSIONS: FrozenSet[str] = frozenset({".consensusxml"})


def _normalized_suffix(path: PathLike) -> str:
    """Return the relevant suffix, handling gzip-compressed files."""

    suffixes = [suffix.lower() for suffix in Path(path).suffixes]
    if not suffixes:
        return ""

    if suffixes[-1] == ".gz" and len(suffixes) >= 2:
        return suffixes[-2]

    return suffixes[-1]


def ensure_allowed_suffix(path: PathLike, allowed: Iterable[str], kind: str) -> str:
    """Validate the suffix for a given kind of OpenMS object."""

    normalized = _normalized_suffix(path)
    allowed_set = {suffix.lower() for suffix in allowed}

    if normalized not in allowed_set:
        allowed_str = ", ".join(sorted(allowed_set)) or "<none>"
        raise ValueError(
            f"Unsupported {kind} file extension for '{path}'. "
            f"Supported extensions: {allowed_str}"
        )

    return normalized

