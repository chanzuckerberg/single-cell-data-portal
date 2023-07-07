import pytest

from backend.curation.api.v1.curation.collections.common import _with_duplicates_removed


@pytest.fixture
def sample_asset_list():
    return [
        {"filetype": "csv", "name": "data.csv"},
        {"filetype": "json", "name": "data.json"},
        {"filetype": "csv", "name": "data2.csv"},
        {"filetype": "txt", "name": "data.txt"},
        {"filetype": "json", "name": "data2.json"},
    ]


def test_remove_duplicate_dataset_assets(sample_asset_list):
    expected_result = [
        {"filetype": "csv", "name": "data.csv"},
        {"filetype": "json", "name": "data.json"},
        {"filetype": "txt", "name": "data.txt"},
    ]
    result = _with_duplicates_removed(sample_asset_list)
    assert result == expected_result


@pytest.fixture
def sample_asset_list_without_duplicates():
    return [
        {"filetype": "csv", "name": "data.csv"},
        {"filetype": "json", "name": "data.json"},
        {"filetype": "txt", "name": "data.txt"},
    ]


def test_remove_duplicate_dataset_assets_no_duplicates(sample_asset_list_without_duplicates):
    expected_result = [
        {"filetype": "csv", "name": "data.csv"},
        {"filetype": "json", "name": "data.json"},
        {"filetype": "txt", "name": "data.txt"},
    ]
    result = _with_duplicates_removed(sample_asset_list_without_duplicates)
    assert result == expected_result
