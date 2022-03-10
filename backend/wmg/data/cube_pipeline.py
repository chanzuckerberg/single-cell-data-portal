import logging
import sys
from typing import List

from backend.wmg.data.extract import copy_datasets_to_instance, get_s3_uris
from backend.wmg.data.load import update_s3_resources
from backend.wmg.data.transform import get_cells_by_tissue_type, generate_cell_ordering
from backend.wmg.data.wmg_cube import create_cube

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def load(dataset_directory: List, group_name: str, validate: bool):
    """
    Given the path to a directory containing one or more h5ad files and a group name, call the h5ad loading function
    on all files, loading/concatenating the datasets together under the group name
    """
    pass


def load_data_and_create_cube(path_to_datasets: str, corpus_name: str, log_level):
    """
    Function to copy H5AD datasets (from a preconfiugred s3 bucket) to the path given then,
    open, transform, normalize and concatenate them together as a tiledb object with a global gene index
    under the given corpus name.
    A cube of summary statistics across genes (but queryable by the given dimensions/data) is then generated
    and stored under big-cube.
    ## TODO add function to get cell count totals
    A per tissue mapping of cell ontologies is generated and the files are copied to s3 under a shared timestamp,.
    On success the least recent set of files are removed from s3
    """
    s3_uris = get_s3_uris()
    copy_datasets_to_instance(s3_uris, path_to_datasets)
    load(path_to_datasets, corpus_name)
    create_cube(corpus_name)  # Todo dry up with create cube func in fixtures
    cell_type_by_tissue = get_cells_by_tissue_type(corpus_name)
    generate_cell_ordering(cell_type_by_tissue)
    update_s3_resources()
    return True


if __name__ == "__main__":
    rv = load_data_and_create_cube()
    sys.exit(rv)
