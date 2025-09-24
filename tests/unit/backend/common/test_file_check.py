from backend.layers.common.files_check import check_file


def test_valid_file(tmp_path):
    file = tmp_path / "valid.txt"
    file.write_bytes(b"data")
    assert check_file(str(file)) is True


def test_empty_file(tmp_path):
    file = tmp_path / "empty.txt"
    file.write_bytes(b"")
    assert check_file(str(file)) is False


def test_directory(tmp_path):
    assert check_file(str(tmp_path)) is False


def test_nonexistent_file(tmp_path):
    missing = tmp_path / "missing.txt"
    assert check_file(str(missing)) is False
