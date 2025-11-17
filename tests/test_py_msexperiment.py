import numpy as np
import pandas as pd
import pytest

oms = pytest.importorskip("pyopenms")

from openms_python.io import read_mzml, write_mzml, stream_mzml
from openms_python.py_msexperiment import Py_MSExperiment
import openms_python.py_msexperiment as py_ms_module
from openms_python.py_msspectrum import Py_MSSpectrum


def build_experiment(num_spectra: int = 3) -> Py_MSExperiment:
    experiment = oms.MSExperiment()
    for idx in range(num_spectra):
        spec = oms.MSSpectrum()
        spec.setRT(10.0 + idx)
        ms_level = 1 if idx < num_spectra - 1 else 2
        spec.setMSLevel(ms_level)
        spec.setNativeID(f"scan={idx}")
        mz = np.array([100.0 + idx, 150.0 + idx, 200.0 + idx])
        intensity = np.array([10.0 * (idx + 1), 20.0 * (idx + 1), 5.0 * (idx + 1)])
        spec.set_peaks((mz.tolist(), intensity.tolist()))
        if ms_level == 2:
            precursor = oms.Precursor()
            precursor.setMZ(400.0)
            precursor.setCharge(2)
            spec.setPrecursors([precursor])
        experiment.addSpectrum(spec)
    return Py_MSExperiment(experiment)


def test_py_msexperiment_iteration_and_summary():
    exp = build_experiment()
    assert len(exp) == 3
    assert exp.rt_range == (10.0, 12.0)
    assert exp.ms_levels == {1, 2}

    ms1 = list(exp.ms1_spectra())
    ms2 = list(exp.ms2_spectra())
    assert len(ms1) == 2
    assert len(ms2) == 1
    assert ms2[0].precursor_mz == pytest.approx(400.0)

    summary = exp.summary()
    assert summary["nr_spectra"] == 3
    assert summary["n_ms1_spectra"] == 2
    assert summary["n_ms2_spectra"] == 1
    assert summary["total_peaks"] == 9


def test_py_msexperiment_pythonic_indexing():
    exp = build_experiment(5)

    assert exp[0].native_id == "scan=0"
    assert exp[-1].native_id == "scan=4"

    every_other = exp[::2]
    assert isinstance(every_other, Py_MSExperiment)
    assert len(every_other) == 3
    assert [spec.native_id for spec in every_other] == ["scan=0", "scan=2", "scan=4"]

    reversed_exp = exp[::-1]
    assert [spec.native_id for spec in reversed_exp] == [f"scan={idx}" for idx in range(4, -1, -1)]

    empty_slice = exp[10:20]
    assert len(empty_slice) == 0


def test_py_msexperiment_dataframe_conversion_and_filters():
    exp = build_experiment()

    spectra_df = exp.to_dataframe(include_peaks=False)
    assert set(["spectrum_index", "native_id", "retention_time", "ms_level", "n_peaks"]).issubset(
        spectra_df.columns
    )
    assert spectra_df.shape[0] == 3

    peaks_df = exp.to_dataframe(include_peaks=True, ms_level=1)
    assert peaks_df["ms_level"].unique().tolist() == [1]
    assert peaks_df.shape[0] == 6  # two spectra * three peaks

    mz_filtered = exp.filter_by_mz_range(100.0, 160.0)
    assert all(len(spec) == 2 for spec in mz_filtered)

    rt_filtered = exp.filter_by_rt_range(10.0, 11.0)
    assert len(rt_filtered) == 2

    sliced = exp[0:2]
    assert len(sliced) == 2
    assert sliced[0].retention_time == pytest.approx(10.0)

    reversed_idx = exp[-1]
    assert isinstance(reversed_idx, Py_MSSpectrum)
    assert reversed_idx.ms_level == 2


def test_py_msexperiment_construction_from_dataframe():
    df = pd.DataFrame(
        {
            "spectrum_id": [0, 0, 1, 1],
            "mz": [100.0, 150.0, 200.0, 250.0],
            "intensity": [10.0, 20.0, 30.0, 40.0],
            "retention_time": [5.0, 5.0, 15.0, 15.0],
            "ms_level": [1, 1, 2, 2],
            "native_id": ["scan=0", "scan=0", "scan=1", "scan=1"],
        }
    )

    exp = Py_MSExperiment.from_dataframe(df)
    assert len(exp) == 2
    assert exp[0].retention_time == pytest.approx(5.0)
    assert exp[1].ms_level == 2


def test_mzml_roundtrip(tmp_path):
    exp = build_experiment()
    output = tmp_path / "roundtrip.mzML"

    write_mzml(exp, output)
    assert output.exists()

    loaded = read_mzml(output)
    assert len(loaded) == len(exp)
    assert [spec.ms_level for spec in loaded] == [spec.ms_level for spec in exp]
    assert [spec.native_id for spec in loaded] == [spec.native_id for spec in exp]
    assert loaded[0].retention_time == pytest.approx(exp[0].retention_time)


def test_stream_mzml_iterates_spectra(tmp_path):
    exp = build_experiment(4)
    path = tmp_path / "streaming.mzML"

    write_mzml(exp, path)

    with stream_mzml(path) as spectra:
        native_ids = [spec.native_id for spec in spectra]

    assert native_ids == [f"scan={idx}" for idx in range(4)]


def test_stream_mzml_native_spectra(tmp_path):
    exp = build_experiment(2)
    path = tmp_path / "stream-native.mzML"

    write_mzml(exp, path)

    with stream_mzml(path, as_wrapper=False) as spectra:
        first = next(spectra)
        assert isinstance(first, oms.MSSpectrum)

    # Ensure context manager drained without iterating all spectra
    assert path.exists()


def test_py_msexperiment_load_store_by_extension(tmp_path):
    exp = build_experiment()
    output = tmp_path / "extension-test.mzML"

    exp.store(output)
    assert output.exists()

    loaded = Py_MSExperiment().load(output)
    assert len(loaded) == len(exp)
    assert [spec.native_id for spec in loaded] == [spec.native_id for spec in exp]


def test_py_msexperiment_rejects_unknown_extension(tmp_path):
    exp = build_experiment()
    bad_path = tmp_path / "experiment.txt"

    with pytest.raises(ValueError):
        exp.store(bad_path)

    with pytest.raises(ValueError):
        Py_MSExperiment().load(bad_path)


def test_peak_picking_convenience(monkeypatch):
    exp = build_experiment()
    original_ms1 = exp[0].mz.tolist()
    processed = []

    class DummyPicker:
        def __init__(self):
            self._params = oms.Param()

        def getParameters(self):
            return self._params

        def setParameters(self, params):
            self._params = params

        def pick(self, source, dest):
            processed.append(source.getNativeID())
            dest.set_peaks(([source.getRT()], [source.size()]))

    monkeypatch.setitem(py_ms_module.PEAK_PICKER_REGISTRY, "hires", DummyPicker)

    picked = exp.pick_peaks(params={"signal_to_noise": 5.0})

    assert processed == ["scan=0", "scan=1"]
    assert picked[0].mz.shape == (1,)
    assert picked[0].mz[0] == pytest.approx(exp[0].retention_time)
    assert picked[2].mz.tolist() == exp[2].mz.tolist()
    assert exp[0].mz.tolist() == original_ms1


def test_peak_picking_inplace_all_levels(monkeypatch):
    exp = build_experiment()
    processed = []

    class DummyPicker:
        def __init__(self):
            self._params = oms.Param()

        def getParameters(self):
            return self._params

        def setParameters(self, params):
            self._params = params

        def pick(self, source, dest):
            processed.append(source.getNativeID())
            dest.set_peaks(([42.0], [84.0]))

    monkeypatch.setitem(py_ms_module.PEAK_PICKER_REGISTRY, "hires", DummyPicker)

    result = exp.pick_peaks(ms_levels=None, inplace=True)

    assert result is exp
    assert len(processed) == len(exp)
    assert all(np.allclose(spec.mz, 42.0) for spec in exp)

