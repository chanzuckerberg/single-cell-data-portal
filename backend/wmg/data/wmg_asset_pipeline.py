import logging
import pathlib
import sys
import time
import tiledb

from backend.atlas_asset_pipelines.integrated_corpus import extract
from backend.atlas_asset_pipelines.integrated_corpus.job import build_integrated_corpus
from backend.atlas_asset_pipelines.expression_summary_cube.job import create_expression_summary_cube
from backend.atlas_asset_pipelines.cell_count_cube.job import create_cell_count_cube
from backend.atlas_asset_pipelines.ontology_ordering.job import generate_cell_ordering
from backend.wmg.data.transform import generate_primary_filter_dimensions
from backend.wmg.data.update_snapshot import update_s3_resources
from backend.wmg.data.schemas.corpus_schema import create_tdb

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def wmg_asset_pipeline(
    path_to_h5ad_datasets: str, corpus_name: str = "corpus_group", snapshot_path=None, extract_data=True
):
    """
    Function to copy H5AD datasets (from a preconfiugred s3 bucket) to the path given then,
    open, transform, normalize and concatenate them together as a tiledb object with a global gene index
    under the given integrated_corpus name.
    A cube of expression summary statistics (known as the expression summary cube) across genes (but queryable by
    the given dimensions/data) is then generated and stored under big-cube.
    A per-tissue mapping of cell ontologies is generated and the files are copied to s3 under a shared timestamp.
    On success the least recent set of files are removed from s3
    """
    if not snapshot_path:
        timestamp = int(time.time())
        snapshot_path = f"{pathlib.Path().resolve()}/{timestamp}"
    corpus_path = f"{snapshot_path}/{corpus_name}"

    if extract_data:
        s3_uris = extract.get_dataset_s3_uris()
        extract.copy_datasets_to_instance(s3_uris, path_to_h5ad_datasets)
        logger.info("Copied datasets to instance")

    if not tiledb.VFS().is_dir(corpus_path):
        create_tdb(snapshot_path, corpus_name)

    build_integrated_corpus(dataset_directory=path_to_h5ad_datasets, corpus_path=corpus_path, validate=True)
    logger.info("Loaded datasets into integrated_corpus")

    create_expression_summary_cube(corpus_path)
    logger.info("Built expression summary cube")

    generate_primary_filter_dimensions(snapshot_path, corpus_name, timestamp)
    create_cell_count_cube(corpus_path)
    logger.info("Built cell count cube")

    generate_cell_ordering(snapshot_path, corpus_path)
    logger.info("Generated cell ordering json file")

    update_s3_resources(snapshot_path, timestamp)
    logger.info("Copied snapshot to s3")

    # todo invalidate distribution id in cloudfront


if __name__ == "__main__":
    wmg_asset_pipeline("datasets", ".")
    sys.exit()
