import shutil
import unittest
from unittest.mock import Mock, patch

import pandas as pd
import tiledb

from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME, EXPRESSION_SUMMARY_CUBE_NAME
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
)
from backend.wmg.pipeline.expression_summary_and_cell_counts import create_expression_summary_and_cell_counts_cubes
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import sort_dataframe
from tests.test_utils.mocks import MockCensusParameters, mock_return_dataset_dict_w_publications
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class ExpressionSummaryAndCellCountsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with (
            tiledb.open(f"{cls.temp_cube_dir.name}/{EXPRESSION_SUMMARY_CUBE_NAME}") as es_cube,
            tiledb.open(f"{cls.temp_cube_dir.name}/{CELL_COUNTS_CUBE_NAME}") as cc_cube,
        ):
            cls.expected_expression_summary_df = sort_dataframe(es_cube.df[:])
            cls.expected_cell_counts_df = sort_dataframe(cc_cube.df[:])
        shutil.rmtree(f"{cls.temp_cube_dir.name}/{EXPRESSION_SUMMARY_CUBE_NAME}")
        shutil.rmtree(f"{cls.temp_cube_dir.name}/{CELL_COUNTS_CUBE_NAME}")

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

        shutil.rmtree(f"{self.temp_cube_dir.name}/{EXPRESSION_SUMMARY_CUBE_NAME}")
        shutil.rmtree(f"{self.temp_cube_dir.name}/{CELL_COUNTS_CUBE_NAME}")

    def test_expression_summary_and_cell_counts(self):
        with (
            patch(
                "backend.wmg.pipeline.expression_summary_and_cell_counts.CensusParameters",
                new=MockCensusParameters,
            ),
            patch("backend.wmg.pipeline.expression_summary.tiledb.consolidate", new=Mock()),
            patch("backend.wmg.pipeline.expression_summary.tiledb.vacuum", new=Mock()),
            patch(
                "backend.wmg.pipeline.cell_counts.return_dataset_dict_w_publications",
                new=mock_return_dataset_dict_w_publications,
            ),
            patch(
                "backend.wmg.pipeline.expression_summary.return_dataset_dict_w_publications",
                new=mock_return_dataset_dict_w_publications,
            ),
        ):
            create_expression_summary_and_cell_counts_cubes(self.temp_cube_dir.name)

        with (
            tiledb.open(f"{self.temp_cube_dir.name}/{EXPRESSION_SUMMARY_CUBE_NAME}") as es_cube,
            tiledb.open(f"{self.temp_cube_dir.name}/{CELL_COUNTS_CUBE_NAME}") as cc_cube,
        ):
            expression_summary_df = sort_dataframe(es_cube.df[:])
            cell_counts_df = sort_dataframe(cc_cube.df[:])
            pd.testing.assert_frame_equal(expression_summary_df, self.expected_expression_summary_df)
            pd.testing.assert_frame_equal(cell_counts_df, self.expected_cell_counts_df)

        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG))
