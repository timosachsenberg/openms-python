"""Pythonic wrapper for pyOpenMS FeatureMap objects."""
from __future__ import annotations

from typing import Iterator, Optional, Union

import pyopenms as oms


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
            index = key
            if index < 0:
                index += len(self)
            if index < 0 or index >= len(self):
                raise IndexError(f"Feature index {key} out of range [0, {len(self)})")
            return self._feature_map[index]

        raise TypeError(f"Invalid index type: {type(key)}")

    def append(self, feature: oms.Feature) -> None:
        """Append a :class:`pyopenms.Feature` to the map."""

        self._feature_map.push_back(feature)
