"""Pythonic wrapper for pyOpenMS FeatureMap objects."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, TYPE_CHECKING, Union

import numpy as np
import pyopenms as oms

from ._io_utils import ensure_allowed_suffix, FEATURE_MAP_EXTENSIONS
from .py_feature import Py_Feature

if TYPE_CHECKING:  # pragma: no cover - typing helper
    from .py_identifications import Identifications


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

    def append(self, feature: Any) -> 'Py_FeatureMap':
        """Append a :class:`pyopenms.Feature` to the map."""

        self._feature_map.push_back(self._as_native_feature(feature))
        return self

    def extend(self, features: Iterable[Any]) -> 'Py_FeatureMap':
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

    # ==================== High-level Workflows ====================

    def annotate_with_identifications(
        self,
        identifications: 'Identifications',
        *,
        mz_tolerance: float = 20.0,
        rt_tolerance: float = 5.0,
        ppm: bool = True,
    ) -> 'Py_FeatureMap':
        """Return a copy annotated with identifications within tolerance."""

        from .py_identifications import Identifications

        if not isinstance(identifications, Identifications):
            raise TypeError("identifications must be an Identifications instance")

        annotated = oms.FeatureMap()

        for feature in self._feature_map:
            assigned: List[oms.PeptideIdentification] = []
            for peptide in identifications.peptides.native:
                if peptide.hasRT() and abs(float(peptide.getRT()) - feature.getRT()) > rt_tolerance:
                    continue
                if peptide.hasMZ() and not self._within_mz(feature.getMZ(), float(peptide.getMZ()), mz_tolerance, ppm):
                    continue
                assigned.append(oms.PeptideIdentification(peptide))

            new_feature = oms.Feature(feature)
            new_feature.setPeptideIdentifications(assigned)
            annotated.push_back(new_feature)

        return Py_FeatureMap(annotated)

    @classmethod
    def align(
        cls,
        feature_maps: Sequence[Union['Py_FeatureMap', oms.FeatureMap]],
        *,
        reference_index: int = 0,
    ) -> List['Py_FeatureMap']:
        """Shift RT coordinates so that all maps share the same median RT."""

        if not feature_maps:
            return []

        py_maps = [cls._ensure_py_feature_map(fmap) for fmap in feature_maps]
        if reference_index < 0 or reference_index >= len(py_maps):
            raise IndexError("reference_index out of range")

        def median_rt(fmap: oms.FeatureMap) -> float:
            rts = [feature.getRT() for feature in fmap]
            return float(np.median(rts)) if rts else 0.0

        reference_median = median_rt(py_maps[reference_index].native)

        aligned: List[Py_FeatureMap] = []
        for fmap in py_maps:
            shift = reference_median - median_rt(fmap.native)
            adjusted = oms.FeatureMap()
            for feature in fmap.native:
                shifted = oms.Feature(feature)
                shifted.setRT(feature.getRT() + shift)
                adjusted.push_back(shifted)
            aligned.append(Py_FeatureMap(adjusted))

        return aligned

    @classmethod
    def link(
        cls,
        feature_maps: Sequence[Union['Py_FeatureMap', oms.FeatureMap]],
        *,
        mz_tolerance: float = 0.01,
        rt_tolerance: float = 5.0,
    ) -> 'Py_ConsensusMap':
        """Link feature maps into a simplified :class:`Py_ConsensusMap`."""

        from .py_consensusmap import Py_ConsensusMap

        if not feature_maps:
            return Py_ConsensusMap()

        py_maps = [cls._ensure_py_feature_map(fmap) for fmap in feature_maps]
        clusters: List[Dict[str, Union[float, List[Tuple[int, oms.Feature]]]]] = []

        for map_index, fmap in enumerate(py_maps):
            for feature in fmap.native:
                match = None
                for cluster in clusters:
                    if (
                        abs(cluster["mz"] - feature.getMZ()) <= mz_tolerance
                        and abs(cluster["rt"] - feature.getRT()) <= rt_tolerance
                    ):
                        match = cluster
                        break
                if match is None:
                    match = {"mz": feature.getMZ(), "rt": feature.getRT(), "members": []}
                    clusters.append(match)
                match["members"].append((map_index, feature))
                mz_values = [member[1].getMZ() for member in match["members"]]
                rt_values = [member[1].getRT() for member in match["members"]]
                match["mz"] = float(np.mean(mz_values))
                match["rt"] = float(np.mean(rt_values))

        consensus = oms.ConsensusMap()
        headers: Dict[int, oms.ColumnHeader] = {}
        for idx, fmap in enumerate(py_maps):
            header = oms.ColumnHeader()
            header.label = f"map_{idx}"
            header.size = fmap.native.size()
            headers[idx] = header
        consensus.setColumnHeaders(headers)

        for cluster in clusters:
            cf = oms.ConsensusFeature()
            cf.setMZ(float(cluster["mz"]))
            cf.setRT(float(cluster["rt"]))
            intensities = []
            for map_index, feature in cluster["members"]:
                peak = oms.Peak2D()
                peak.setRT(feature.getRT())
                peak.setMZ(feature.getMZ())
                peak.setIntensity(feature.getIntensity())
                cf.insert(map_index, peak, 0)
                intensities.append(feature.getIntensity())
            if intensities:
                cf.setIntensity(float(np.sum(intensities)))
            consensus.push_back(cf)

        return Py_ConsensusMap(consensus)

    # ==================== Private Helpers ====================

    def _as_native_feature(self, feature: Any) -> oms.Feature:
        """Return a native :class:`pyopenms.Feature` from supported inputs."""

        if isinstance(feature, oms.Feature):
            return feature

        native = getattr(feature, "native", None)
        if isinstance(native, oms.Feature):
            return native

        raise TypeError(
            "feature entries must be pyopenms.Feature or expose a native Feature via `.native`"
        )

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
    def _within_mz(feature_mz: float, id_mz: float, tolerance: float, ppm: bool) -> bool:
        delta = abs(feature_mz - id_mz)
        allowed = tolerance * feature_mz / 1e6 if ppm else tolerance
        return delta <= allowed

    @staticmethod
    def _ensure_py_feature_map(feature_map: Union['Py_FeatureMap', oms.FeatureMap]) -> 'Py_FeatureMap':
        if isinstance(feature_map, Py_FeatureMap):
            return feature_map
        if isinstance(feature_map, oms.FeatureMap):
            return Py_FeatureMap(feature_map)
        raise TypeError("feature_maps must contain Py_FeatureMap or pyopenms.FeatureMap")
    def _as_native_feature(feature: Union[oms.Feature, Py_Feature]) -> oms.Feature:
        if isinstance(feature, Py_Feature):
            return feature.native
        if isinstance(feature, oms.Feature):
            return feature
        raise TypeError(
            "Features must be pyopenms.Feature or Py_Feature instances, "
            f"got {type(feature).__name__}"
        )
