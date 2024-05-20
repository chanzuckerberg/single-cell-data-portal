import shutil
import unittest

import pandas as pd
import tiledb

from backend.wmg.data.snapshot import (
    CELL_COUNTS_DIFFEXP_CUBE_NAME,
    EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME,
    EXPRESSION_SUMMARY_DIFFEXP_SIMPLE_CUBE_NAME,
)
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_DIFFEXP_CUBES_CREATED_FLAG,
)
from backend.wmg.pipeline.expression_summary_and_cell_counts_diffexp import (
    create_expression_summary_and_cell_counts_diffexp_cubes,
)
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import sort_dataframe
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class ExpressionSummaryDiffexpTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with tiledb.open(f"{cls.temp_cube_dir.name}/{EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME}") as es_cube:
            setattr(cls, f"expected_{EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME}_df", sort_dataframe(es_cube.df[:]))
        with tiledb.open(f"{cls.temp_cube_dir.name}/{CELL_COUNTS_DIFFEXP_CUBE_NAME}") as cell_counts_cube:
            setattr(cls, f"expected_{CELL_COUNTS_DIFFEXP_CUBE_NAME}_df", sort_dataframe(cell_counts_cube.df[:]))
        with tiledb.open(f"{cls.temp_cube_dir.name}/{EXPRESSION_SUMMARY_DIFFEXP_SIMPLE_CUBE_NAME}") as es_simple_cube:
            setattr(
                cls, f"expected_{EXPRESSION_SUMMARY_DIFFEXP_SIMPLE_CUBE_NAME}_df", sort_dataframe(es_simple_cube.df[:])
            )

        for cube_name in [
            EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME,
            EXPRESSION_SUMMARY_DIFFEXP_SIMPLE_CUBE_NAME,
            CELL_COUNTS_DIFFEXP_CUBE_NAME,
        ]:
            shutil.rmtree(f"{cls.temp_cube_dir.name}/{cube_name}")

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def setUp(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG] = True
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(EXPRESSION_SUMMARY_AND_CELL_COUNTS_DIFFEXP_CUBES_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)
        for cube_name in [
            EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME,
            EXPRESSION_SUMMARY_DIFFEXP_SIMPLE_CUBE_NAME,
            CELL_COUNTS_DIFFEXP_CUBE_NAME,
        ]:
            shutil.rmtree(f"{self.temp_cube_dir.name}/{cube_name}")

    def test_expression_summary_diffexp(self):
        create_expression_summary_and_cell_counts_diffexp_cubes(self.temp_cube_dir.name)
        for cube_name in [
            EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME,
            EXPRESSION_SUMMARY_DIFFEXP_SIMPLE_CUBE_NAME,
            CELL_COUNTS_DIFFEXP_CUBE_NAME,
        ]:
            with tiledb.open(f"{self.temp_cube_dir.name}/{cube_name}") as es_cube:
                expression_summary_diffexp_df = sort_dataframe(es_cube.df[:])
                pd.testing.assert_frame_equal(expression_summary_diffexp_df, getattr(self, f"expected_{cube_name}_df"))

        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_DIFFEXP_CUBES_CREATED_FLAG))
