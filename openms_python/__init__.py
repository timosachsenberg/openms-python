"""
openms_python: A Pythonic wrapper around pyOpenMS

This package provides a more intuitive, Python-friendly interface to OpenMS
for mass spectrometry data analysis.

Example:
    >>> from openms_python import MSExperiment
    >>> exp = MSExperiment.from_file('data.mzML')
    >>> print(f"Loaded {len(exp)} spectra")
    >>> for spec in exp.ms1_spectra():
    ...     print(f"RT: {spec.retention_time:.2f}, Peaks: {len(spec)}")
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Union

__version__ = "0.1.2"
__author__ = "MiniMax Agent"

from .py_msexperiment import Py_MSExperiment
from .py_msspectrum import Py_MSSpectrum
from .py_feature import Py_Feature
from .py_featuremap import Py_FeatureMap
from .py_consensusmap import Py_ConsensusMap
from .py_identifications import (
    ProteinIdentifications,
    PeptideIdentifications,
    Identifications,
)
from .workflows import (
    DigestedPeptide,
    TheoreticalSpectrumRecord,
    ProteinStream,
    PeptideStream,
    SpectrumStream,
    stream_theoretical_spectra_from_fasta,
    map_identifications_to_features,
    align_feature_maps,
    link_features,
    export_quant_table,
)
from .io import read_mzml, write_mzml, stream_mzml

_EXAMPLE_PACKAGE = "openms_python.examples"


def get_example(name: str, *, load: bool = False, target_dir: Union[str, Path, None] = None) -> Union[str, bytes]:
    """Return a bundled example file path or its contents.

    Parameters
    ----------
    name:
        File name inside ``openms_python/examples`` (e.g. ``"small.mzML"``).
    load:
        When ``True`` the file is returned as raw bytes instead of a file path.
    target_dir:
        Optional directory where the example file should be copied before returning the
        path. This is useful when the package is installed as a zipped archive because
        :mod:`importlib.resources` may otherwise expose a temporary path.

    Returns
    -------
    str or bytes
        The absolute path to the requested example file or its binary content when
        ``load`` is ``True``.
    """

    if load:
        try:
            return resources.read_binary(_EXAMPLE_PACKAGE, name)
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"Example file '{name}' not found") from exc

    if target_dir is not None:
        destination = Path(target_dir).expanduser().resolve()
        destination.mkdir(parents=True, exist_ok=True)
        destination_file = destination / name
        try:
            data = resources.read_binary(_EXAMPLE_PACKAGE, name)
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"Example file '{name}' not found") from exc
        destination_file.write_bytes(data)
        return str(destination_file)

    try:
        with resources.path(_EXAMPLE_PACKAGE, name) as resource_path:
            resolved = Path(resource_path)
            if not resolved.exists():
                raise FileNotFoundError(name)
            return str(resolved)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Example file '{name}' not found") from exc

__all__ = [
    "Py_MSExperiment",
    "Py_MSSpectrum",
    "Py_Feature",
    "Py_FeatureMap",
    "Py_ConsensusMap",
    "ProteinIdentifications",
    "PeptideIdentifications",
    "Identifications",
    "DigestedPeptide",
    "TheoreticalSpectrumRecord",
    "ProteinStream",
    "PeptideStream",
    "SpectrumStream",
    "stream_theoretical_spectra_from_fasta",
    "map_identifications_to_features",
    "align_feature_maps",
    "link_features",
    "export_quant_table",
    "read_mzml",
    "write_mzml",
    "stream_mzml",
    "get_example",
]
