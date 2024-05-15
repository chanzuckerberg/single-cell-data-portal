import json
import unittest
from unittest.mock import patch

from backend.cellguide.pipeline.explorer_cxgs import get_valid_cxgs
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils import compare_dicts
from tests.test_utils.mocks import mock_get_folders_from_s3
from tests.unit.backend.cellguide.pipeline.constants import (
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    VALID_EXPLORER_CXGS_FIXTURE_FILENAME,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class TestValidExplorerCxgs(unittest.TestCase):
    def test__valid_explorer_cxgs(self):
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{VALID_EXPLORER_CXGS_FIXTURE_FILENAME}", "r") as f:
            expected__valid_cxgs = json.load(f)
        with patch(
            "backend.cellguide.pipeline.explorer_cxgs.get_folders_from_s3",
            new=mock_get_folders_from_s3,
        ):
            valid_cxgs = get_valid_cxgs()
        self.assertTrue(compare_dicts(convert_dataclass_to_dict_and_strip_nones(valid_cxgs), expected__valid_cxgs))
