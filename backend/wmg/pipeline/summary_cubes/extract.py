import numpy as np
import pandas as pd
import tiledb

from backend.wmg.data.schemas.corpus_schema import OBS_ARRAY_NAME, VAR_ARRAY_NAME


def extract_var_data(tdb_group: str, ctx: tiledb.Ctx) -> list:
    """
    Extract var (gene) data from the concatenated corpus
    """
    with tiledb.open(f"{tdb_group}/{VAR_ARRAY_NAME}", ctx=ctx) as var:
        gene_ontology_term_ids = var.query(dims=["gene_ontology_term_id"], attrs=["var_idx"], use_arrow=False).df[:]
        gene_ontology_term_ids.sort_values(by="var_idx", inplace=True)
        return gene_ontology_term_ids


def extract_obs_data(tdb_group: str, cube_dims: list) -> pd.DataFrame:
    """
    Extract obs (cell) data from the concatenated corpus
    """
    with tiledb.open(f"{tdb_group}/{OBS_ARRAY_NAME}") as obs:
        obs_df = obs.query(use_arrow=False).df[:]

        # filter out observations in the 'filter_cells' attribute
        cell_labels = obs_df[np.logical_not(obs_df["filter_cells"])]

    cell_labels.sort_values(by=["obs_idx"], inplace=True, ignore_index=True)

    cell_labels = pd.DataFrame(
        data={k: cell_labels[k].astype("category") for k in cube_dims},
        index=cell_labels.index,
    )
    return cell_labels
