import json
import unittest

from backend.wmg.data.snapshot import PRIMARY_FILTER_DIMENSIONS_FILENAME
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
    PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG,
)
from backend.wmg.pipeline.primary_filter_dimensions import create_primary_filter_dimensions
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class PrimaryFilterDimensionsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with open(f"{cls.temp_cube_dir.name}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}") as f:
            cls.expected_primary_filter_dimensions = json.load(f)

    def setUp(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG] = True
        pipeline_state[EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG] = True
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def test_primary_filter_dimensions(self):
        create_primary_filter_dimensions(self.temp_cube_dir.name)
        with open(f"{self.temp_cube_dir.name}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}") as f:
            primary_filter_dimensions = json.load(f)
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG))
        self.assertTrue(compare_dicts(primary_filter_dimensions, self.expected_primary_filter_dimensions))
