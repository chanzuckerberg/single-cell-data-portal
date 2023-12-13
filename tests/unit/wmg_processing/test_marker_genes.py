import shutil
import unittest
from unittest.mock import patch

import pandas as pd
import tiledb

from backend.wmg.data.snapshot import MARKER_GENES_CUBE_NAME
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
    MARKER_GENES_CUBE_CREATED_FLAG,
    PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG,
)
from backend.wmg.pipeline.marker_genes import create_marker_genes_cube
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import sort_dataframe
from tests.test_utils.mocks import mock_bootstrap_rows_percentiles
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class MarkerGenesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with tiledb.open(f"{cls.temp_cube_dir.name}/{MARKER_GENES_CUBE_NAME}") as mg_cube:
            cls.expected_marker_gene_df = sort_dataframe(mg_cube.df[:])
        shutil.rmtree(f"{cls.temp_cube_dir.name}/{MARKER_GENES_CUBE_NAME}")

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def setUp(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG] = True
        pipeline_state[PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG] = True
        pipeline_state[EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG] = True
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(MARKER_GENES_CUBE_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)
        shutil.rmtree(f"{self.temp_cube_dir.name}/{MARKER_GENES_CUBE_NAME}")

    def test_marker_genes_cube(self):
        with patch(
            "backend.cellguide.pipeline.computational_marker_genes.computational_markers.bootstrap_rows_percentiles",
            new=mock_bootstrap_rows_percentiles,
        ):
            create_marker_genes_cube(self.temp_cube_dir.name)

        with tiledb.open(f"{self.temp_cube_dir.name}/{MARKER_GENES_CUBE_NAME}") as mg_cube:
            marker_genes_df = sort_dataframe(mg_cube.df[:])
            pd.testing.assert_frame_equal(marker_genes_df, self.expected_marker_gene_df)

        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(MARKER_GENES_CUBE_CREATED_FLAG))
