import argparse
import os
import shutil
import sys
import tempfile
from enum import Enum
from unittest.mock import Mock, patch

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
from tests.test_utils.mocks import (
    MockCensusParameters,
    mock_bootstrap_rows_percentiles,
    mock_get_datasets_from_curation_endpoint,
    mock_return_dataset_dict_w_publications,
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


class FixtureType(Enum):
    all = "all"
    cell_counts = "cell_counts"
    expression_summary = "expression_summary"
    marker_genes = "marker_genes"
    expression_summary_default = "expression_summary_default"
    dataset_metadata = "dataset_metadata"
    primary_filter_dimensions = "primary_filter_dimensions"
    cell_type_orderings = "cell_type_orderings"
    filter_relationships = "filter_relationships"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "new_snapshot",
        help="New snapshot location",
    )
    parser.add_argument(
        "--fixture_type",
        help="Type of fixture to generate",
        choices=[member.value for member in FixtureType],
        default="all",
    )
    args = parser.parse_args()

    new_snapshot = args.new_snapshot
    fixture_type = args.fixture_type

    if not os.path.isdir(new_snapshot):
        os.mkdir(new_snapshot)

    with (
        tempfile.TemporaryDirectory() as temp_dir,
        patch(
            "backend.wmg.pipeline.expression_summary_and_cell_counts.CensusParameters",
            new=MockCensusParameters,
        ),
        patch(
            "backend.cellguide.pipeline.computational_marker_genes.computational_markers.bootstrap_rows_percentiles",
            new=mock_bootstrap_rows_percentiles,
        ),
        patch(
            "backend.wmg.pipeline.dataset_metadata.get_datasets_from_discover_api",
            new=mock_get_datasets_from_curation_endpoint,
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
        corpus_path = os.path.join(temp_dir, os.path.basename(os.path.normpath(new_snapshot)))
        run_pipeline(corpus_path=corpus_path, skip_validation=True)

        for cube_name in [
            CELL_COUNTS_CUBE_NAME,
            MARKER_GENES_CUBE_NAME,
            EXPRESSION_SUMMARY_CUBE_NAME,
            EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME,
        ]:
            if cube_name == fixture_type or fixture_type == FixtureType.all.value:
                with tiledb.open(os.path.join(corpus_path, cube_name)) as cube:
                    df = cube.df[:]
                path = os.path.join(new_snapshot, cube_name + ".csv")
                print(f"Writing {path}")
                df.to_csv(path)
                os.system(f"rm -rf {path}.gz")
                os.system(f"gzip {path}")

        for filename in [
            DATASET_METADATA_FILENAME,
            PRIMARY_FILTER_DIMENSIONS_FILENAME,
            CELL_TYPE_ORDERINGS_FILENAME,
            FILTER_RELATIONSHIPS_FILENAME,
        ]:
            if filename.split(".json")[0] == fixture_type or fixture_type == FixtureType.all.value:
                path = os.path.join(new_snapshot, filename)
                print(f"Writing {path}")
                shutil.copy(os.path.join(corpus_path, filename), path)
                os.system(f"rm -rf {path}.gz")
                os.system(f"gzip {path}")
