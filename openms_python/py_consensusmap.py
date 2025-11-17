"""Pythonic wrapper for pyOpenMS ConsensusMap objects."""
from __future__ import annotations

from pathlib import Path
from typing import Iterator, Optional, Union

import pyopenms as oms

from ._io_utils import ensure_allowed_suffix, CONSENSUS_MAP_EXTENSIONS


class Py_ConsensusMap:
    """Provide sequence-like access and IO helpers for ConsensusMap objects."""

    def __init__(self, native_map: Optional[oms.ConsensusMap] = None):
        self._consensus_map = native_map if native_map is not None else oms.ConsensusMap()

    @property
    def native(self) -> oms.ConsensusMap:
        """Return the underlying :class:`pyopenms.ConsensusMap`."""

        return self._consensus_map

    def __len__(self) -> int:  # pragma: no cover - trivial
        return int(self._consensus_map.size())

    def __iter__(self) -> Iterator[oms.ConsensusFeature]:
        for index in range(len(self)):
            yield self[index]

    def __getitem__(self, key: Union[int, slice]) -> Union[oms.ConsensusFeature, 'Py_ConsensusMap']:
        """Return consensus features or a sliced :class:`Py_ConsensusMap`."""

        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            sliced_map = oms.ConsensusMap()
            for idx in range(start, stop, step):
                sliced_map.push_back(self._consensus_map[idx])
            return Py_ConsensusMap(sliced_map)

        if isinstance(key, int):
            index = key
            if index < 0:
                index += len(self)
            if index < 0 or index >= len(self):
                raise IndexError(f"Consensus feature index {key} out of range [0, {len(self)})")
            return self._consensus_map[index]

        raise TypeError(f"Invalid index type: {type(key)}")

    def load(self, filepath: Union[str, Path]) -> 'Py_ConsensusMap':
        """Load a consensus map from disk by inspecting the extension."""

        ensure_allowed_suffix(filepath, CONSENSUS_MAP_EXTENSIONS, "ConsensusMap")
        oms.ConsensusXMLFile().load(str(filepath), self._consensus_map)
        return self

    def store(self, filepath: Union[str, Path]) -> 'Py_ConsensusMap':
        """Store the consensus map to disk, validating the output extension."""

        ensure_allowed_suffix(filepath, CONSENSUS_MAP_EXTENSIONS, "ConsensusMap")
        oms.ConsensusXMLFile().store(str(filepath), self._consensus_map)
        return self
