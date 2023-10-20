import unittest
from unittest.mock import patch

import tiledb

from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME, EXPRESSION_SUMMARY_CUBE_NAME
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    WMG_PIPELINE_TEST_RUN_KEY,
)
from backend.wmg.pipeline.expression_summary_and_cell_counts import create_expression_summary_and_cell_counts_cubes
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import TempEnvironmentVariable
from tests.test_utils.mocks import MockExpressionSummaryBuilder, mock_create_cell_counts_cube
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class ExpressionSummaryAndCellCountsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with (
            tiledb.open(f"{cls.temp_cube_dir.name}/{EXPRESSION_SUMMARY_CUBE_NAME}") as es_cube,
            tiledb.open(f"{cls.temp_cube_dir.name}/{CELL_COUNTS_CUBE_NAME}") as cc_cube,
        ):
            cls.expected_expression_summary_df = es_cube.df[:]
            cls.cell_counts_df = cc_cube.df[:]
        # shutil.rmtree(f"{cls.temp_cube_dir.name}/{EXPRESSION_SUMMARY_CUBE_NAME}")
        # shutil.rmtree(f"{cls.temp_cube_dir.name}/{CELL_COUNTS_CUBE_NAME}")

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

        # shutil.rmtree(f"{self.temp_cube_dir.name}/{EXPRESSION_SUMMARY_CUBE_NAME}")
        # shutil.rmtree(f"{self.temp_cube_dir.name}/{CELL_COUNTS_CUBE_NAME}")

    def test_expression_summary_and_cell_counts(self):
        with (
            TempEnvironmentVariable(WMG_PIPELINE_TEST_RUN_KEY, "true"),
            patch(
                "backend.wmg.pipeline.expression_summary_and_cell_counts.ExpressionSummaryCubeBuilder",
                new=MockExpressionSummaryBuilder,
            ),
            patch(
                "backend.wmg.pipeline.expression_summary_and_cell_counts.create_cell_counts_cube",
                new=mock_create_cell_counts_cube,
            ),
        ):
            create_expression_summary_and_cell_counts_cubes(self.temp_cube_dir.name)

        # TODO: figure out how to mock cellxgene_census
        # with (
        #     tiledb.open(f"{self.temp_cube_dir.name}/{EXPRESSION_SUMMARY_CUBE_NAME}") as es_cube,
        #     tiledb.open(f"{self.temp_cube_dir.name}/{CELL_COUNTS_CUBE_NAME}") as cc_cube,
        # ):
        #     expression_summary_df = es_cube.df[:]
        #     cell_counts_df = cc_cube.df[:]
        #     pd.testing.assert_frame_equal(expression_summary_df, self.expected_expression_summary_df)
        #     pd.testing.assert_frame_equal(cell_counts_df, self.cell_counts_df)

        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG))
