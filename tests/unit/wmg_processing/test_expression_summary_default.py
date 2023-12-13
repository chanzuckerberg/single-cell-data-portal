import shutil
import unittest

import pandas as pd
import tiledb

from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
)
from backend.wmg.pipeline.expression_summary_default import create_expression_summary_default_cube
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import sort_dataframe
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class ExpressionSummaryDefaultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with (tiledb.open(f"{cls.temp_cube_dir.name}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}") as es_cube,):
            cls.expected_expression_summary_default_df = sort_dataframe(es_cube.df[:])
        shutil.rmtree(f"{cls.temp_cube_dir.name}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}")

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def setUp(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG] = True
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)
        shutil.rmtree(f"{self.temp_cube_dir.name}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}")

    def test_expression_summary_default(self):
        create_expression_summary_default_cube(self.temp_cube_dir.name)

        with (tiledb.open(f"{self.temp_cube_dir.name}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}") as es_cube,):
            expression_summary_default_df = sort_dataframe(es_cube.df[:])
            pd.testing.assert_frame_equal(expression_summary_default_df, self.expected_expression_summary_default_df)

        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG))
