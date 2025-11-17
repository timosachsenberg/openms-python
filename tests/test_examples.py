from pathlib import Path

import pytest

from openms_python import get_example


def test_get_example_returns_existing_path():
    example_path = get_example("small.mzML")
    path = Path(example_path)
    assert path.exists()
    assert path.name == "small.mzML"
    assert path.read_text().startswith("<?xml")


def test_get_example_load_and_target_dir(tmp_path):
    example_bytes = get_example("small.mzML", load=True)
    assert example_bytes.startswith(b"<?xml")

    destination = get_example("small.mzML", target_dir=tmp_path)
    dest_path = Path(destination)
    assert dest_path.parent == tmp_path
    assert dest_path.read_bytes() == example_bytes


def test_get_example_missing_file():
    with pytest.raises(FileNotFoundError):
        get_example("does-not-exist.mzML")
