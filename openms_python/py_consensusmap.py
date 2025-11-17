"""Pythonic wrapper for pyOpenMS ConsensusMap objects."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, Optional, Union

import pandas as pd
import pyopenms as oms

from ._io_utils import ensure_allowed_suffix, CONSENSUS_MAP_EXTENSIONS


class Py_ConsensusMap:
    """Provide sequence-like access and IO helpers for ConsensusMap objects."""

    _BASE_DF_COLUMNS = [
        "unique_id",
        "rt",
        "mz",
        "intensity",
        "charge",
        "width",
        "quality",
    ]

    _OPTIONAL_FIELD_SETTERS = {
        "unique_id": ("setUniqueId", int),
        "charge": ("setCharge", int),
        "width": ("setWidth", float),
        "quality": ("setQuality", float),
    }

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
            index = self._normalize_index(key)
            return self._consensus_map[index]

        raise TypeError(f"Invalid index type: {type(key)}")

    def load(self, filepath: Union[str, Path]) -> 'Py_ConsensusMap':
        """Load a consensus map from disk by inspecting the extension."""

        ensure_allowed_suffix(filepath, CONSENSUS_MAP_EXTENSIONS, "ConsensusMap")
        oms.ConsensusXMLFile().load(str(filepath), self._consensus_map)
        return self

    def append(self, feature: oms.ConsensusFeature) -> 'Py_ConsensusMap':
        """Append a :class:`pyopenms.ConsensusFeature` to the map."""

        self._consensus_map.push_back(feature)
        return self

    def extend(self, features: Iterable[oms.ConsensusFeature]) -> 'Py_ConsensusMap':
        """Append multiple consensus features to the map."""

        for feature in features:
            self.append(feature)
        return self

    def remove(self, index: int) -> 'Py_ConsensusMap':
        """Remove the consensus feature at *index* and return ``self``."""

        normalized = self._normalize_index(index)
        self._delete_indices([normalized])
        return self

    def __delitem__(self, key: Union[int, slice]) -> None:
        """Delete one or more consensus features using Python's semantics."""

        if isinstance(key, int):
            self.remove(key)
            return

        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            indices = list(range(start, stop, step))
            self._delete_indices(indices)
            return

        raise TypeError(f"Invalid deletion index type: {type(key)}")

    def store(self, filepath: Union[str, Path]) -> 'Py_ConsensusMap':
        """Store the consensus map to disk, validating the output extension."""

        ensure_allowed_suffix(filepath, CONSENSUS_MAP_EXTENSIONS, "ConsensusMap")
        oms.ConsensusXMLFile().store(str(filepath), self._consensus_map)
        return self

    # ==================== pandas integration ====================

    def to_dataframe(self) -> pd.DataFrame:
        """Convert consensus features to a :class:`pandas.DataFrame`."""

        rows = []
        for feature in self._consensus_map:
            entry = {
                "unique_id": feature.getUniqueId(),
                "rt": feature.getRT(),
                "mz": feature.getMZ(),
                "intensity": feature.getIntensity(),
                "charge": feature.getCharge(),
                "width": feature.getWidth(),
                "quality": feature.getQuality(),
            }

            meta_keys: list = []
            feature.getKeys(meta_keys)
            for raw_key in meta_keys:
                key = raw_key.decode() if isinstance(raw_key, (bytes, bytearray)) else str(raw_key)
                entry[key] = feature.getMetaValue(key)

            rows.append(entry)

        if not rows:
            return pd.DataFrame(columns=self._BASE_DF_COLUMNS)

        df = pd.DataFrame(rows)
        ordered = self._BASE_DF_COLUMNS + [
            column for column in df.columns if column not in self._BASE_DF_COLUMNS
        ]
        return df[ordered]

    def get_df(self) -> pd.DataFrame:
        """Compatibility alias for :meth:`to_dataframe`."""

        return self.to_dataframe()

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> 'Py_ConsensusMap':
        """Construct a :class:`Py_ConsensusMap` from a :class:`pandas.DataFrame`."""

        required = {"rt", "mz", "intensity"}
        missing = required - set(df.columns)
        if missing:
            missing_str = ", ".join(sorted(missing))
            raise ValueError(f"DataFrame is missing required columns: {missing_str}")

        consensus_map = oms.ConsensusMap()
        meta_columns = [
            column
            for column in df.columns
            if column not in required and column not in cls._OPTIONAL_FIELD_SETTERS
        ]

        for record in df.to_dict(orient="records"):
            feature = oms.ConsensusFeature()
            feature.setRT(float(record["rt"]))
            feature.setMZ(float(record["mz"]))
            feature.setIntensity(float(record["intensity"]))

            for column, (setter_name, caster) in cls._OPTIONAL_FIELD_SETTERS.items():
                if column not in record:
                    continue
                value = record[column]
                if cls._is_missing(value):
                    continue
                setter = getattr(feature, setter_name)
                setter(cls._cast_value(value, caster))

            for column in meta_columns:
                value = record.get(column)
                if column in required or cls._is_missing(value):
                    continue
                feature.setMetaValue(column, cls._normalize_meta_value(value))

            consensus_map.push_back(feature)

        return cls(consensus_map)

    @classmethod
    def from_df(cls, df: pd.DataFrame) -> 'Py_ConsensusMap':
        """Alias for :meth:`from_dataframe` matching :meth:`get_df`."""

        return cls.from_dataframe(df)

    # ==================== Private Helpers ====================

    def _normalize_index(self, index: int) -> int:
        length = len(self)
        if length == 0:
            raise IndexError("ConsensusMap is empty")
        if index < 0:
            index += length
        if index < 0 or index >= length:
            raise IndexError(f"Consensus feature index {index} out of range [0, {length})")
        return index

    def _delete_indices(self, indices: Iterable[int]) -> None:
        drop = sorted(set(indices))
        if not drop:
            return

        length = len(self)
        source_map = self._consensus_map
        drop_set = set(drop)

        new_map = oms.ConsensusMap(source_map)
        new_map.clear(False)

        for idx in range(length):
            if idx in drop_set:
                continue
            new_map.push_back(oms.ConsensusFeature(source_map[idx]))

        self._consensus_map = new_map

    @staticmethod
    def _is_missing(value: object) -> bool:
        try:
            return pd.isna(value)
        except TypeError:
            return False

    @staticmethod
    def _cast_value(value: object, caster):
        if hasattr(value, "item"):
            try:
                value = value.item()
            except Exception:
                pass
        return caster(value)

    @staticmethod
    def _normalize_meta_value(value: object):
        if hasattr(value, "item"):
            try:
                return value.item()
            except Exception:
                return value
        return value
