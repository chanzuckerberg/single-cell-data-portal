import json

from tests.unit.cellguide_pipeline.constants import (
    ASCTB_MASTER_SHEET_FIXTURE_FILENAME,
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
)


def mock_get_asctb_master_sheet():
    with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ASCTB_MASTER_SHEET_FIXTURE_FILENAME}", "r") as f:
        return json.load(f)
