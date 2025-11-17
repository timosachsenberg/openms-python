import pytest

oms = pytest.importorskip("pyopenms")

from openms_python.py_feature import Py_Feature
from openms_python.py_msspectrum import Py_MSSpectrum


def test_py_feature_meta_behaves_like_dict():
    feature = Py_Feature()

    feature["label"] = "Sample1"
    assert feature["label"] == "Sample1"
    assert feature.get("label") == "Sample1"
    assert "label" in feature

    feature.setdefault("condition", "control")
    assert feature["condition"] == "control"
    feature.setdefault("condition", "changed")
    assert feature["condition"] == "control"

    feature.update({"score": 42}, replicate="r1")
    assert feature["score"] == 42
    assert feature["replicate"] == "r1"

    popped = feature.pop("score")
    assert popped == 42
    with pytest.raises(KeyError):
        _ = feature["score"]

    feature.clear()
    assert feature.get("label") is None


def test_py_msspectrum_meta_access():
    spectrum = Py_MSSpectrum(oms.MSSpectrum())

    spectrum["IonInjectTime"] = 13.5
    assert spectrum["IonInjectTime"] == pytest.approx(13.5)
    assert spectrum.get("missing", 0.0) == 0.0

    spectrum.update([(b"Instrument", "QExactive")])
    assert spectrum["Instrument"] == "QExactive"

    value = spectrum.pop("Instrument", None)
    assert value == "QExactive"
    assert spectrum.get("Instrument") is None
