import logging
import os

import tiledb
from tiledbsoma import ExperimentAxisQuery

from backend.wmg.data.schemas.cube_schema import (
    cell_counts_logical_dims,
    cell_counts_schema,
)
from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
)
from backend.wmg.data.utils import log_func_runtime
from backend.wmg.pipeline.constants import (
    DIMENSION_NAME_MAP_CENSUS_TO_WMG,
)
from backend.wmg.pipeline.utils import (
    create_empty_cube_if_needed,
    remove_accents,
    return_dataset_dict_w_publications,
)

logger = logging.getLogger(__name__)


@log_func_runtime
def create_cell_counts_cube(*, query: ExperimentAxisQuery, corpus_path: str):
    """
    Create cell count cube and write to disk
    """
    obs_df = query.obs().concat().to_pandas()
    obs_df = obs_df.rename(columns=DIMENSION_NAME_MAP_CENSUS_TO_WMG)

    logger.info("Creating the cell counts cube and filter relationships graph.")

    df = (
        obs_df.groupby(
            by=[dim for dim in cell_counts_logical_dims if dim != "publication_citation"],
            as_index=False,
        ).size()
    ).rename(columns={"size": "n_cells"})

    dataset_dict = return_dataset_dict_w_publications()
    df["publication_citation"] = [
        remove_accents(dataset_dict.get(dataset_id, "No Publication")) for dataset_id in df["dataset_id"]
    ]
    n_cells = df["n_cells"].to_numpy()
    df["n_cells"] = n_cells

    uri = os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)
    create_empty_cube_if_needed(uri, cell_counts_schema)
    logger.info("Writing cell counts cube.")
    tiledb.from_pandas(uri, df, mode="append")
