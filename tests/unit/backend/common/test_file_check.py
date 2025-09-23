import os
import tempfile

from backend.layers.common.files_check import check_file


def test_valid_file():
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp.write(b"data")
        tmp_path = tmp.name
    assert check_file(tmp_path) is True
    os.remove(tmp_path)


def test_empty_file():
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp_path = tmp.name
    assert check_file(tmp_path) is False
    os.remove(tmp_path)


def test_unreadable_file():
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp.write(b"data")
        tmp_path = tmp.name
    os.chmod(tmp_path, 0o000)
    assert check_file(tmp_path) is False
    os.chmod(tmp_path, 0o644)
    os.remove(tmp_path)


def test_directory():
    with tempfile.TemporaryDirectory() as tmpdir:
        assert check_file(tmpdir) is False


def test_nonexistent_file():
    assert check_file("/tmp/nonexistent_file_12345") is False
