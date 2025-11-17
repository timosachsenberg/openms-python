import numpy as np
import pandas as pd
import pytest

oms = pytest.importorskip("pyopenms")

from openms_python.py_msspectrum import Py_MSSpectrum


def create_native_spectrum():
    spec = oms.MSSpectrum()
    spec.setRT(12.3)
    spec.setMSLevel(2)
    spec.setNativeID("controllerType=0 controllerNumber=1 scan=15")
    spec.set_peaks((
        [100.5, 150.2, 300.4, 320.5],
        [200.0, 500.0, 50.0, 1000.0],
    ))
    precursor = oms.Precursor()
    precursor.setMZ(543.21)
    precursor.setCharge(3)
    spec.setPrecursors([precursor])
    return spec


def test_py_msspectrum_properties_and_peak_operations():
    wrapper = Py_MSSpectrum(create_native_spectrum())

    assert pytest.approx(wrapper.retention_time) == 12.3
    wrapper.retention_time = 14.6
    assert pytest.approx(wrapper.retention_time) == 14.6

    assert wrapper.ms_level == 2
    wrapper.ms_level = 1
    assert wrapper.is_ms1
    assert not wrapper.is_ms2
    wrapper.ms_level = 2  # restore original

    assert wrapper.precursor_mz == pytest.approx(543.21)
    assert wrapper.precursor_charge == 3
    assert wrapper.scan_number == 15
    assert len(wrapper) == 4

    mz, intensity = wrapper.peaks
    assert np.allclose(mz, np.array([100.5, 150.2, 300.4, 320.5]))
    assert np.allclose(intensity, np.array([200.0, 500.0, 50.0, 1000.0]))

    assert wrapper.total_ion_current == pytest.approx(1750.0)
    assert wrapper.base_peak_mz == pytest.approx(320.5)
    assert wrapper.base_peak_intensity == pytest.approx(1000.0)

    filtered = wrapper.filter_by_mz(120.0, 310.0)
    assert np.all(filtered.mz >= 120.0) and np.all(filtered.mz <= 310.0)
    assert len(filtered) == 2

    filtered_intensity = wrapper.filter_by_intensity(400.0)
    assert np.all(filtered_intensity.intensity >= 400.0)
    assert len(filtered_intensity) == 2

    top_one = wrapper.top_n_peaks(1)
    assert len(top_one) == 1
    assert top_one.base_peak_mz == pytest.approx(320.5)

    normalized = wrapper.normalize_intensity(max_value=100.0)
    assert normalized.base_peak_intensity == pytest.approx(100.0)

    new_mz = np.array([50.0, 75.0])
    new_intensity = np.array([10.0, 20.0])
    wrapper.peaks = (new_mz, new_intensity)
    assert np.allclose(wrapper.mz, new_mz)
    assert np.allclose(wrapper.intensity, new_intensity)


def test_py_msspectrum_dataframe_helpers_round_trip():
    df = pd.DataFrame({"mz": [100.0, 200.0], "intensity": [10.0, 20.0]})

    wrapper = Py_MSSpectrum.from_dataframe(df, retention_time=5.5, ms_level=1, native_id="scan=1")

    recreated_df = wrapper.to_dataframe()
    assert recreated_df.shape == (2, 2)
    assert recreated_df["mz"].tolist() == [100.0, 200.0]
    assert recreated_df["intensity"].tolist() == [10.0, 20.0]
    assert wrapper.retention_time == pytest.approx(5.5)
    assert wrapper.ms_level == 1
    assert wrapper.native_id == "scan=1"

