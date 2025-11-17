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

__version__ = "0.1.2"
__author__ = "MiniMax Agent"

from .py_msexperiment import Py_MSExperiment
from .py_msspectrum import Py_MSSpectrum
from .py_featuremap import Py_FeatureMap
from .io import read_mzml, write_mzml

__all__ = [
    "Py_MSExperiment",
    "Py_MSSpectrum",
    "Py_FeatureMap",
    "read_mzml",
    "write_mzml",
]
