import os
import shutil
import sys
import tempfile
from unittest.mock import patch

import tiledb

from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    CELL_TYPE_ORDERINGS_FILENAME,
    DATASET_METADATA_FILENAME,
    EXPRESSION_SUMMARY_CUBE_NAME,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME,
    FILTER_RELATIONSHIPS_FILENAME,
    MARKER_GENES_CUBE_NAME,
    PRIMARY_FILTER_DIMENSIONS_FILENAME,
)
from backend.wmg.pipeline import run_pipeline
from backend.wmg.pipeline.constants import WMG_PIPELINE_TEST_RUN_KEY
from tests.test_utils import TempEnvironmentVariable
from tests.test_utils.mocks import (
    mock_bootstrap_rows_percentiles,
    mock_get_datasets_from_curation_endpoint,
)

# Add the root directory to the Python module search path so you can reference backend
# without needing to move this script to the root directory to run it.
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)

"""
This generates realistic test data by invoking the pipeline on a tiny subset of cells.

You can run this script from the root directory of single-cell-data-portal:
python -m scripts.generate_realistic_test_snapshot tests/unit/backend/wmg/fixtures/realistic-test-snapshot
"""

if __name__ == "__main__":
    new_snapshot = sys.argv[1]

    if not os.path.isdir(new_snapshot):
        os.mkdir(new_snapshot)

    with (
        TempEnvironmentVariable(WMG_PIPELINE_TEST_RUN_KEY, "true"),
        tempfile.TemporaryDirectory() as temp_dir,
        patch(
            "backend.cellguide.pipeline.computational_marker_genes.computational_markers.bootstrap_rows_percentiles",
            new=mock_bootstrap_rows_percentiles,
        ),
        patch(
            "backend.wmg.pipeline.dataset_metadata.get_datasets_from_discover_api",
            new=mock_get_datasets_from_curation_endpoint,
        ),
    ):
        corpus_path = os.path.join(temp_dir, "test-snapshot")
        run_pipeline(corpus_path=corpus_path)

        for cube_name in [
            CELL_COUNTS_CUBE_NAME,
            MARKER_GENES_CUBE_NAME,
            EXPRESSION_SUMMARY_CUBE_NAME,
            EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME,
        ]:
            with tiledb.open(os.path.join(corpus_path, cube_name)) as cube:
                df = cube.df[:]
            df.to_csv(os.path.join(new_snapshot, cube_name + ".csv"))

        for filename in [
            DATASET_METADATA_FILENAME,
            PRIMARY_FILTER_DIMENSIONS_FILENAME,
            CELL_TYPE_ORDERINGS_FILENAME,
            FILTER_RELATIONSHIPS_FILENAME,
        ]:
            shutil.copy(os.path.join(corpus_path, filename), os.path.join(new_snapshot, filename))

        os.system(f"gzip -r {new_snapshot}/*")
