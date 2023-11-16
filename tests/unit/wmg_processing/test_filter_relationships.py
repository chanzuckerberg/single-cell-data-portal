import json
import unittest

from backend.wmg.data.snapshot import FILTER_RELATIONSHIPS_FILENAME
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    FILTER_RELATIONSHIPS_CREATED_FLAG,
)
from backend.wmg.pipeline.filter_relationships import create_filter_relationships_graph
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import CompareDictsAddin
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class FilterRelationshipsTests(unittest.TestCase, CompareDictsAddin):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with open(f"{cls.temp_cube_dir.name}/{FILTER_RELATIONSHIPS_FILENAME}") as f:
            cls.expected_filter_relationships = json.load(f)

    def setUp(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG] = True
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(FILTER_RELATIONSHIPS_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def test_filter_relationships(self):
        create_filter_relationships_graph(self.temp_cube_dir.name)
        with open(f"{self.temp_cube_dir.name}/{FILTER_RELATIONSHIPS_FILENAME}") as f:
            filter_relationships = json.load(f)
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(FILTER_RELATIONSHIPS_CREATED_FLAG))
        self.assert_dicts_equal(filter_relationships, self.expected_filter_relationships)
