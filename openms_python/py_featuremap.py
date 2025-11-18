"""Pythonic wrapper for pyOpenMS FeatureMap objects."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Union

import pyopenms as oms
import pandas as pd
from ._io_utils import ensure_allowed_suffix, FEATURE_MAP_EXTENSIONS
from .py_feature import Py_Feature
from .py_identifications import Identifications


class Py_FeatureMap:
    """Provide sequence-like access to :class:`pyopenms.FeatureMap`."""

    _BASE_DF_COLUMNS = [
        "unique_id",
        "rt",
        "mz",
        "intensity",
        "charge",
        "width",
        "overall_quality",
    ]

    _OPTIONAL_FIELD_SETTERS = {
        "unique_id": ("setUniqueId", int),
        "charge": ("setCharge", int),
        "width": ("setWidth", float),
        "overall_quality": ("setOverallQuality", float),
    }

    def __init__(self, native_map: Optional[oms.FeatureMap] = None):
        self._feature_map = native_map if native_map is not None else oms.FeatureMap()

    @property
    def native(self) -> oms.FeatureMap:
        """Return the underlying :class:`pyopenms.FeatureMap`."""

        return self._feature_map

    def __len__(self) -> int:  # pragma: no cover - trivial
        return int(self._feature_map.size())

    def __iter__(self) -> Iterator[Py_Feature]:
        for index in range(len(self)):
            yield self[index]

    def __getitem__(self, key: Union[int, slice]) -> Union[Py_Feature, 'Py_FeatureMap']:
        """Return individual features or a sliced :class:`Py_FeatureMap`."""

        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            sliced_map = oms.FeatureMap()
            for idx in range(start, stop, step):
                sliced_map.push_back(self._feature_map[idx])
            return Py_FeatureMap(sliced_map)

        if isinstance(key, int):
            index = self._normalize_index(key)
            return Py_Feature(self._feature_map[index])

        raise TypeError(f"Invalid index type: {type(key)}")

    def append(self, feature: Union[oms.Feature, Py_Feature]) -> 'Py_FeatureMap':
        """Append a :class:`pyopenms.Feature` to the map."""

        self._feature_map.push_back(self._as_native_feature(feature))
        return self

    def extend(self, features: Iterable[Union[oms.Feature, Py_Feature]]) -> 'Py_FeatureMap':
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

    # ==================== Protein inference ====================

    def infer_proteins(
        self,
        *,
        algorithm: str = "basic",
        params: Optional[Union[oms.Param, Dict[str, Union[int, float, str]]]] = None,
        include_unassigned: bool = False,
        greedy_group_resolution: bool = True,
        experimental_design: Optional[oms.ExperimentalDesign] = None,
    ) -> Identifications:
        """Run protein inference on the identifications stored in this map."""

        identifications = self._collect_identifications()
        return identifications.infer_proteins(
            algorithm=algorithm,
            params=params,
            include_unassigned=include_unassigned,
            greedy_group_resolution=greedy_group_resolution,
            experimental_design=experimental_design,
        )

    # ==================== pandas integration ====================

    def to_dataframe(self) -> pd.DataFrame:
        """Convert all features into a :class:`pandas.DataFrame`."""

        rows = []
        for feature in self._feature_map:
            entry = {
                "unique_id": feature.getUniqueId(),
                "rt": feature.getRT(),
                "mz": feature.getMZ(),
                "intensity": feature.getIntensity(),
                "charge": feature.getCharge(),
                "width": feature.getWidth(),
                "overall_quality": feature.getOverallQuality(),
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
        ordered_columns = self._BASE_DF_COLUMNS + [
            column for column in df.columns if column not in self._BASE_DF_COLUMNS
        ]
        return df[ordered_columns]

    def get_df(self) -> pd.DataFrame:
        """Compatibility alias for :meth:`to_dataframe`."""

        return self.to_dataframe()

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> 'Py_FeatureMap':
        """Create a :class:`Py_FeatureMap` from a :class:`pandas.DataFrame`."""

        required = {"rt", "mz", "intensity"}
        missing = required - set(df.columns)
        if missing:
            missing_str = ", ".join(sorted(missing))
            raise ValueError(f"DataFrame is missing required columns: {missing_str}")

        feature_map = oms.FeatureMap()
        meta_columns = [
            column for column in df.columns if column not in required and column not in cls._OPTIONAL_FIELD_SETTERS
        ]

        for record in df.to_dict(orient="records"):
            feature = oms.Feature()
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

            feature_map.push_back(feature)

        return cls(feature_map)

    @classmethod
    def from_df(cls, df: pd.DataFrame) -> 'Py_FeatureMap':
        """Alias for :meth:`from_dataframe` matching :meth:`get_df`."""

        return cls.from_dataframe(df)

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

    def _as_native_feature(self, feature: Union[oms.Feature, Py_Feature]) -> oms.Feature:
        if isinstance(feature, Py_Feature):
            return feature.native
        if isinstance(feature, oms.Feature):
            return feature
        raise TypeError(
            "Features must be pyopenms.Feature or Py_Feature instances, "
            f"got {type(feature).__name__}"
        )

    def _collect_identifications(self) -> Identifications:
        proteins = self._feature_map.getProteinIdentifications()
        peptides = list(self._feature_map.getUnassignedPeptideIdentifications())
        peptides.extend(self._feature_map.get_assigned_peptide_identifications())
        return Identifications(proteins, peptides)
