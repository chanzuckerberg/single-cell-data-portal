import logging
import pathlib
import sys
import time
import tiledb

from backend.corpus_asset_pipelines import integrated_corpus
from backend.corpus_asset_pipelines import expression_summary_cube

from backend.wmg.data.load_cube import upload_artifacts_to_s3, make_snapshot_active
from backend.wmg.data.schemas.corpus_schema import create_tdb
from backend.wmg.data.transform import (
    generate_primary_filter_dimensions,
    get_cell_types_by_tissue,
    generate_cell_ordering,
)
from backend.wmg.data.wmg_cube import create_cell_count_cube

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def load_data_and_create_cube(
    path_to_h5ad_datasets: str,
    corpus_name: str = "corpus_group",
    snapshot_path=None,
    extract_data=True,
    validate_cube=True,
):
    """
    Function to copy H5AD datasets (from a preconfiugred s3 bucket) to the path given then,
    open, transform, normalize and concatenate them together as a tiledb object with a global gene index
    under the given corpus name.
    A cube of expression summary statistics (known as the expression summary cube) across genes (but queryable by
    the given dimensions/data) is then generated and stored under big-cube.
    A per-tissue mapping of cell ontologies is generated and the files are copied to s3 under a shared timestamp,.
    On success the least recent set of files are removed from s3
    """
    if not snapshot_path:
        snapshot_id = int(time.time())
        snapshot_path = f"{pathlib.Path().resolve()}/{snapshot_id}"
    corpus_path = f"{snapshot_path}/{corpus_name}"
    if not tiledb.VFS().is_dir(corpus_path):
        create_tdb(snapshot_path, corpus_name)

    integrated_corpus.run(path_to_h5ad_datasets, corpus_path, extract_data)
    logger.info("Loaded datasets into corpus")
    expression_summary_cube.run(corpus_path, validate_cube)
    logger.info("Build expression summary cube")
    create_cell_count_cube(corpus_path)
    logger.info("Built cell count cube")

    cell_type_by_tissue = get_cell_types_by_tissue(corpus_path)
    generate_cell_ordering(snapshot_path, cell_type_by_tissue)
    generate_primary_filter_dimensions(snapshot_path, corpus_name, snapshot_id)
    logger.info("Generated cell ordering json file")
    upload_artifacts_to_s3(snapshot_path, snapshot_id)
    logger.info("Copied snapshot to s3")
    if validate_cube:
        make_snapshot_active(snapshot_id)
        logger.info(f"Updated latest_snapshot_identifier in s3. Current snapshot id is {snapshot_id}")


if __name__ == "__main__":
    """
    To trigger the batch process in prod via the cli run
    AWS_PROFILE=single-cell-prod aws batch submit-job --job-name $JOB_NAME --job-queue dp-prod --job-definition dp-prod-prodstack-wmg-processing # noqa E501
    """
    # todo pass in validate_cubes as env arg
    load_data_and_create_cube("datasets", ".")
    sys.exit()
