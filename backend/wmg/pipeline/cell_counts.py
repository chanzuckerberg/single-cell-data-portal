import logging
import os

import tiledb
from tiledbsoma import ExperimentAxisQuery

from backend.common.utils.rollup import ancestors
from backend.wmg.data.schemas.cube_schema import (
    cell_counts_logical_dims,
    cell_counts_schema,
)
from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
)
from backend.wmg.pipeline.constants import (
    DIMENSION_NAME_MAP_CENSUS_TO_WMG,
)
from backend.wmg.pipeline.utils import (
    create_empty_cube_if_needed,
    log_func_runtime,
    remove_accents,
    return_dataset_dict_w_publications,
)

logger = logging.getLogger(__name__)


@log_func_runtime
def create_cell_counts_cube(*, query: ExperimentAxisQuery, corpus_path: str, organismId: str):
    """
    Create cell count cube and write to disk
    """
    obs_df = query.obs().concat().to_pandas()
    obs_df = obs_df.rename(columns=DIMENSION_NAME_MAP_CENSUS_TO_WMG)
    obs_df["organism_ontology_term_id"] = organismId

    logger.info("Creating the cell counts cube.")

    df = (
        obs_df.groupby(
            by=[
                dim
                for dim in cell_counts_logical_dims
                if dim != "publication_citation" and dim != "cell_type_ontology_term_id_ancestors"
            ],
            as_index=False,
        ).size()
    ).rename(columns={"size": "n_cells"})

    dataset_dict = return_dataset_dict_w_publications()
    df["publication_citation"] = [
        remove_accents(dataset_dict.get(dataset_id, "No Publication")) for dataset_id in df["dataset_id"]
    ]
    n_cells = df["n_cells"].to_numpy()
    df["n_cells"] = n_cells

    # The following is used to create a dictionary of ancestors for each cell type id.
    # If a cell type id is not already present in the dictionary, it is added along with its ancestors.
    # The 'ancestors' function is used to get the ancestors of a cell type id.
    # The 'apply' function is then used to apply the 'get_ancestors' function to each cell type id in the dataframe.
    # This results in a new column in the dataframe, 'cell_type_ontology_term_id_ancestors', which contains the ancestors of each cell type id.
    # These ancestors are later used as part of the rollup operation.
    logger.info("Creating the cell type ancestors column.")
    ancestors_dict = {}

    def get_ancestors(cell_type_id):
        if cell_type_id not in ancestors_dict:
            ancestors_dict[cell_type_id] = ",".join(ancestors(cell_type_id))
        return ancestors_dict[cell_type_id]

    obs_df["cell_type_ontology_term_id_ancestors"] = obs_df["cell_type_ontology_term_id"].apply(get_ancestors)

    uri = os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)
    create_empty_cube_if_needed(uri, cell_counts_schema)
    logger.info("Writing cell counts cube.")
    tiledb.from_pandas(uri, df, mode="append")
