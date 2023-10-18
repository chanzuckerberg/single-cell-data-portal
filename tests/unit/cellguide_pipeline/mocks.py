import json

from tests.unit.cellguide_pipeline.constants import (
    ASCTB_MASTER_SHEET_FIXTURE_FILENAME,
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    FAKE_COLLECTIONS_FIXTURE_FILENAME,
    FAKE_DATASETS_FIXTURE_FILENAME,
)


def mock_get_asctb_master_sheet():
    with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ASCTB_MASTER_SHEET_FIXTURE_FILENAME}", "r") as f:
        return json.load(f)


def mock_get_datasets_from_curation_endpoint():
    with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{FAKE_DATASETS_FIXTURE_FILENAME}", "r") as f:
        return json.load(f)


def mock_get_collections_from_curation_endpoint():
    with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{FAKE_COLLECTIONS_FIXTURE_FILENAME}", "r") as f:
        return json.load(f)


def mock_get_title_and_citation_from_doi(doi: str):
    return f"title from {doi}\n\n-- citation"
