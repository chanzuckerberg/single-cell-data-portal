import json
import unittest

from backend.wmg.data.snapshot import CELL_TYPE_ANCESTORS_FILENAME
from backend.wmg.pipeline.cell_type_ancestors import create_cell_type_ancestors
from backend.wmg.pipeline.constants import (
    CELL_TYPE_ANCESTORS_CREATED_FLAG,
)
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class CellTypeAncestorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with open(f"{cls.temp_cube_dir.name}/{CELL_TYPE_ANCESTORS_FILENAME}") as f:
            cls.expected_cell_type_ancestors = json.load(f)

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(CELL_TYPE_ANCESTORS_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def test_cell_type_ancestors(self):
        create_cell_type_ancestors(self.temp_cube_dir.name)

        with open(f"{self.temp_cube_dir.name}/{CELL_TYPE_ANCESTORS_FILENAME}") as f:
            cell_type_ancestors = json.load(f)

        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(CELL_TYPE_ANCESTORS_CREATED_FLAG))
        self.assertTrue(compare_dicts(cell_type_ancestors, self.expected_cell_type_ancestors))
