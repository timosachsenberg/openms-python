"""Pythonic wrapper for pyOpenMS FeatureMap objects."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, Optional, Union

import pyopenms as oms
from ._io_utils import ensure_allowed_suffix, FEATURE_MAP_EXTENSIONS


class Py_FeatureMap:
    """Provide sequence-like access to :class:`pyopenms.FeatureMap`."""

    def __init__(self, native_map: Optional[oms.FeatureMap] = None):
        self._feature_map = native_map if native_map is not None else oms.FeatureMap()

    @property
    def native(self) -> oms.FeatureMap:
        """Return the underlying :class:`pyopenms.FeatureMap`."""

        return self._feature_map

    def __len__(self) -> int:  # pragma: no cover - trivial
        return int(self._feature_map.size())

    def __iter__(self) -> Iterator[oms.Feature]:
        for index in range(len(self)):
            yield self[index]

    def __getitem__(self, key: Union[int, slice]) -> Union[oms.Feature, 'Py_FeatureMap']:
        """Return individual features or a sliced :class:`Py_FeatureMap`."""

        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            sliced_map = oms.FeatureMap()
            for idx in range(start, stop, step):
                sliced_map.push_back(self._feature_map[idx])
            return Py_FeatureMap(sliced_map)

        if isinstance(key, int):
            index = self._normalize_index(key)
            return self._feature_map[index]

        raise TypeError(f"Invalid index type: {type(key)}")

    def append(self, feature: oms.Feature) -> 'Py_FeatureMap':
        """Append a :class:`pyopenms.Feature` to the map."""

        self._feature_map.push_back(feature)
        return self

    def extend(self, features: Iterable[oms.Feature]) -> 'Py_FeatureMap':
        """Append multiple features to the map."""

        for feature in features:
            self.append(feature)
        return self

    def remove(self, index: int) -> 'Py_FeatureMap':
        """Remove the feature at *index* and return ``self`` for chaining."""

        normalized = self._normalize_index(index)
        self._delete_indices([normalized])
        return self

    def __delitem__(self, key: Union[int, slice]) -> None:
        """Delete one or multiple features using Python's deletion semantics."""

        if isinstance(key, int):
            self.remove(key)
            return

        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            indices = list(range(start, stop, step))
            self._delete_indices(indices)
            return

        raise TypeError(f"Invalid deletion index type: {type(key)}")

    def load(self, filepath: Union[str, Path]) -> 'Py_FeatureMap':
        """Load a feature map from disk by inspecting the extension."""

        ensure_allowed_suffix(filepath, FEATURE_MAP_EXTENSIONS, "FeatureMap")
        oms.FeatureXMLFile().load(str(filepath), self._feature_map)
        return self

    def store(self, filepath: Union[str, Path]) -> 'Py_FeatureMap':
        """Store the feature map to disk, validating the output extension."""

        ensure_allowed_suffix(filepath, FEATURE_MAP_EXTENSIONS, "FeatureMap")
        oms.FeatureXMLFile().store(str(filepath), self._feature_map)
        return self

    # ==================== Private Helpers ====================

    def _normalize_index(self, index: int) -> int:
        length = len(self)
        if length == 0:
            raise IndexError("FeatureMap is empty")
        if index < 0:
            index += length
        if index < 0 or index >= length:
            raise IndexError(f"Feature index {index} out of range [0, {length})")
        return index

    def _delete_indices(self, indices: Iterable[int]) -> None:
        drop = sorted(set(indices))
        if not drop:
            return

        length = len(self)
        source_map = self._feature_map
        drop_set = set(drop)

        new_map = oms.FeatureMap(source_map)
        new_map.clear(False)

        for idx in range(length):
            if idx in drop_set:
                continue
            new_map.push_back(oms.Feature(source_map[idx]))

        self._feature_map = new_map
