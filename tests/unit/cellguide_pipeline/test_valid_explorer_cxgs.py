import json
import unittest

from backend.cellguide.pipeline.explorer_cxgs import get_valid_cxgs
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils import compare_dicts
from tests.unit.cellguide_pipeline.constants import (
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    VALID_EXPLORER_CXGS_FIXTURE_FILENAME,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class TestMetadataGenerator(unittest.TestCase):
    def test__cell_metadata_generator(self):
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{VALID_EXPLORER_CXGS_FIXTURE_FILENAME}", "r") as f:
            expected__valid_cxgs = json.load(f)

        valid_cxgs = get_valid_cxgs()
        self.assertTrue(compare_dicts(convert_dataclass_to_dict_and_strip_nones(valid_cxgs), expected__valid_cxgs))
