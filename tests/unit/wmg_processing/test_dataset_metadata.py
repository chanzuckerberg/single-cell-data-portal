import json
import unittest
from unittest.mock import patch

from backend.wmg.data.snapshot import DATASET_METADATA_FILENAME
from backend.wmg.pipeline.constants import (
    DATASET_METADATA_CREATED_FLAG,
)
from backend.wmg.pipeline.dataset_metadata import create_dataset_metadata
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import compare_dicts
from tests.test_utils.mocks import mock_get_datasets_from_curation_endpoint
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class DatasetMetadataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with open(f"{cls.temp_cube_dir.name}/{DATASET_METADATA_FILENAME}") as f:
            cls.expected_dataset_metadata = json.load(f)

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(DATASET_METADATA_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def test_dataset_metadata(self):
        with patch(
            "backend.wmg.pipeline.dataset_metadata.get_datasets_from_discover_api",
            new=mock_get_datasets_from_curation_endpoint,
        ):
            create_dataset_metadata(self.temp_cube_dir.name)
        with open(f"{self.temp_cube_dir.name}/{DATASET_METADATA_FILENAME}") as f:
            dataset_metadata = json.load(f)
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(DATASET_METADATA_CREATED_FLAG))
        self.assertTrue(compare_dicts(dataset_metadata, self.expected_dataset_metadata))
