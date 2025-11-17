"""
Pythonic wrapper for pyOpenMS Spectrum classes.
"""

from typing import Tuple, Optional
import numpy as np
import pandas as pd
import pyopenms as oms

from ._meta_mapping import MetaInfoMappingMixin


class Py_MSSpectrum(MetaInfoMappingMixin):
    """
    A Pythonic wrapper around pyOpenMS MSSpectrum.
    
    This class provides intuitive properties and methods for working with
    mass spectra, hiding the verbose C++ API underneath.
    
    Example:
        >>> spec = Spectrum(native_spectrum)
        >>> print(f"RT: {spec.retention_time:.2f} seconds")
        >>> print(f"MS Level: {spec.ms_level}")
        >>> print(f"Number of peaks: {len(spec)}")
        >>> if spec.is_ms1:
        ...     print("This is an MS1 spectrum")
        >>> peaks_df = spec.to_dataframe()
    """
    
    def __init__(self, native_spectrum: oms.MSSpectrum):
        """
        Initialize Spectrum wrapper.
        
        Args:
            native_spectrum: pyOpenMS MSSpectrum object
        """
        self._spectrum = native_spectrum

    # ==================== Meta-info support ====================

    def _meta_object(self) -> oms.MetaInfoInterface:
        return self._spectrum
    
    # ==================== Pythonic Properties ====================
    
    @property
    def retention_time(self) -> float:
        """Get retention time in seconds."""
        return self._spectrum.getRT()
    
    @retention_time.setter
    def retention_time(self, value: float):
        """Set retention time in seconds."""
        self._spectrum.setRT(value)
    
    @property
    def ms_level(self) -> int:
        """Get MS level (1 for MS1, 2 for MS2, etc.)."""
        return self._spectrum.getMSLevel()
    
    @ms_level.setter
    def ms_level(self, value: int):
        """Set MS level."""
        self._spectrum.setMSLevel(value)
    
    @property
    def is_ms1(self) -> bool:
        """Check if this is an MS1 spectrum."""
        return self.ms_level == 1
    
    @property
    def is_ms2(self) -> bool:
        """Check if this is an MS2 spectrum."""
        return self.ms_level == 2
    
    @property
    def precursor_mz(self) -> Optional[float]:
        """Get precursor m/z for MS2+ spectra, None for MS1."""
        if self.ms_level < 2:
            return None
        precursors = self._spectrum.getPrecursors()
        if precursors:
            return precursors[0].getMZ()
        return None
    
    @property
    def precursor_charge(self) -> Optional[int]:
        """Get precursor charge for MS2+ spectra, None for MS1."""
        if self.ms_level < 2:
            return None
        precursors = self._spectrum.getPrecursors()
        if precursors:
            return precursors[0].getCharge()
        return None
    
    @property
    def native_id(self) -> str:
        """Get native ID of the spectrum."""
        return self._spectrum.getNativeID()
    
    @property
    def scan_number(self) -> int:
        """Get scan number (extracted from native ID or -1 if not available)."""
        try:
            # Try to extract scan number from native ID
            native_id = self.native_id
            if "scan=" in native_id:
                return int(native_id.split("scan=")[1].split()[0])
        except (ValueError, IndexError):
            pass
        return -1
    
    @property
    def total_ion_current(self) -> float:
        """Get total ion current (sum of all peak intensities)."""
        _, intensities = self.peaks
        return float(np.sum(intensities))
    
    @property
    def base_peak_mz(self) -> Optional[float]:
        """Get m/z of the base peak (most intense peak)."""
        if len(self) == 0:
            return None
        mz, intensities = self.peaks
        return float(mz[np.argmax(intensities)])
    
    @property
    def base_peak_intensity(self) -> Optional[float]:
        """Get intensity of the base peak."""
        if len(self) == 0:
            return None
        _, intensities = self.peaks
        return float(np.max(intensities))
    
    @property
    def peaks(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get peaks as NumPy arrays.
        
        Returns:
            Tuple of (mz_array, intensity_array)
        """
        mz, intensity = self._spectrum.get_peaks()
        return np.array(mz), np.array(intensity)
    
    @peaks.setter
    def peaks(self, values: Tuple[np.ndarray, np.ndarray]):
        """
        Set peaks from NumPy arrays.
        
        Args:
            values: Tuple of (mz_array, intensity_array)
        """
        mz, intensity = values
        self._spectrum.set_peaks((mz.tolist(), intensity.tolist()))

    @property
    def mz(self) -> np.ndarray:
        """Get m/z values as a NumPy array."""
        mz, _ = self.peaks
        return mz

    @property
    def intensity(self) -> np.ndarray:
        """Get intensity values as a NumPy array."""
        _, intensity = self.peaks
        return intensity
    

    # ==================== Magic Methods ====================
    
    def __len__(self) -> int:
        """Return number of peaks in the spectrum."""
        return self._spectrum.size()
    
    def __repr__(self) -> str:
        """Return string representation."""
        ms_info = f"MS{self.ms_level}"
        if self.precursor_mz is not None:
            ms_info += f" (precursor: {self.precursor_mz:.4f})"
        return (
            f"Spectrum(rt={self.retention_time:.2f}s, {ms_info}, "
            f"peaks={len(self)}, TIC={self.total_ion_current:.2e})"
        )
    
    def __str__(self) -> str:
        """Return human-readable string."""
        return self.__repr__()


    def __iter__(self):
        """Allow dict(self) and list(self) conversions."""
        yield "mz", self.mz.tolist()
        yield "intens", self.intensity.tolist()
        

    
    # ==================== Conversion Methods ====================

    def to_numpy(self) -> np.ndarray:
        """
        Convert spectrum peaks to NumPy arrays.
        
        Returns:
            Tuple of (mz_array, intensity_array)
            
        Example:
            >>> mz, intensity = spec.to_numpy()
        """
        return np.array(self.peaks)

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert spectrum peaks to pandas DataFrame.
        
        Returns:
            DataFrame with columns: mz, intensity
            
        Example:
            >>> df = spec.to_dataframe()
            >>> df.head()
                      mz   intensity
            0  100.0500     1250.5
            1  200.1234     5678.2
            ...
        """
        mz, intensity = self.peaks
        return pd.DataFrame({
            'mz': mz,
            'intensity': intensity
        })
    
    @classmethod
    def from_numpy(mz: np.ndarray, intensity: np.ndarray) -> 'Py_MSSpectrum':
        """
        Create spectrum from NumPy arrays.
        
        Args:
            mz: Array of m/z values
            intensity: Array of intensity values
        """
        spec = oms.MSSpectrum()
        spec.set_peaks((mz.tolist(), intensity.tolist()))
        return Py_MSSpectrum(spec)


    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, **metadata) -> 'Py_MSSpectrum':
        """
        Create spectrum from pandas DataFrame.
        
        Args:
            df: DataFrame with 'mz' and 'intensity' columns
            **metadata: Optional metadata (retention_time, ms_level, etc.)
            
        Returns:
            Spectrum object
            
        Example:
            >>> df = pd.DataFrame({'mz': [100, 200], 'intensity': [50, 100]})
            >>> spec = Spectrum.from_dataframe(df, retention_time=60.5, ms_level=1)
        """
        spec = oms.MSSpectrum()
        spec.set_peaks((df['mz'].values.tolist(), df['intensity'].values.tolist()))
        
        # Set metadata
        if 'retention_time' in metadata:
            spec.setRT(metadata['retention_time'])
        if 'ms_level' in metadata:
            spec.setMSLevel(metadata['ms_level'])
        if 'native_id' in metadata:
            spec.setNativeID(metadata['native_id'])
            
        return cls(spec)
    
    # ==================== Data Manipulation ====================

    def filter_by_mz(self, min_mz: float, max_mz: float) -> 'Py_MSSpectrum':
        """
        Filter peaks by m/z range.
        
        Args:
            min_mz: Minimum m/z value
            max_mz: Maximum m/z value
            
        Returns:
            New Spectrum with filtered peaks
        """
        mz, intensity = self.peaks
        mask = (mz >= min_mz) & (mz <= max_mz)
        
        new_spec = oms.MSSpectrum()
        new_spec.set_peaks((mz[mask].tolist(), intensity[mask].tolist()))
        new_spec.setRT(self.retention_time)
        new_spec.setMSLevel(self.ms_level)
        new_spec.setNativeID(self.native_id)
        
        return Py_MSSpectrum(new_spec)
    
    def filter_by_intensity(self, min_intensity: float) -> 'Py_MSSpectrum':
        """
        Filter peaks by minimum intensity.
        
        Args:
            min_intensity: Minimum intensity threshold
            
        Returns:
            New Spectrum with filtered peaks
        """
        mz, intensity = self.peaks
        mask = intensity >= min_intensity
        
        new_spec = oms.MSSpectrum()
        new_spec.set_peaks((mz[mask].tolist(), intensity[mask].tolist()))
        new_spec.setRT(self.retention_time)
        new_spec.setMSLevel(self.ms_level)
        new_spec.setNativeID(self.native_id)
        
        return Py_MSSpectrum(new_spec)

    def top_n_peaks(self, n: int) -> 'Py_MSSpectrum':
        """
        Keep only the top N most intense peaks.
        
        Args:
            n: Number of peaks to keep
            
        Returns:
            New Spectrum with top N peaks
        """
        mz, intensity = self.peaks
        if len(mz) <= n:
            return self
        
        # Get indices of top N peaks
        top_indices = np.argsort(intensity)[-n:]
        top_indices = np.sort(top_indices)  # Keep m/z order
        
        new_spec = oms.MSSpectrum()
        new_spec.set_peaks((mz[top_indices].tolist(), intensity[top_indices].tolist()))
        new_spec.setRT(self.retention_time)
        new_spec.setMSLevel(self.ms_level)
        new_spec.setNativeID(self.native_id)
        
        return Py_MSSpectrum(new_spec)

    def normalize_intensity(self, max_value: float = 100.0) -> 'Py_MSSpectrum':
        """
        Normalize peak intensities to a maximum value.
        
        Args:
            max_value: Target maximum intensity (default: 100.0)
            
        Returns:
            New Spectrum with normalized intensities
        """
        mz, intensity = self.peaks
        if len(intensity) == 0 or np.max(intensity) == 0:
            return self
        
        normalized = intensity * (max_value / np.max(intensity))
        
        new_spec = oms.MSSpectrum()
        new_spec.set_peaks((mz.tolist(), normalized.tolist()))
        new_spec.setRT(self.retention_time)
        new_spec.setMSLevel(self.ms_level)
        new_spec.setNativeID(self.native_id)

        return Py_MSSpectrum(new_spec)

    # ==================== Access to Native Object ====================
    
    @property
    def native(self) -> oms.MSSpectrum:
        """
        Get the underlying pyOpenMS MSSpectrum object.
        
        Use this when you need to access pyOpenMS-specific methods
        not wrapped by this class.
        """
        return self._spectrum
