import logging

import pandas as pd
import tiledb

from backend.wmg.data.schemas.corpus_schema import OBS_ARRAY_NAME
from backend.wmg.data.schemas.cube_schema import cell_counts_schema
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME
from backend.wmg.data.utils import create_empty_cube

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def extract(corpus_path: str) -> pd.DataFrame:
    return tiledb.open(f"{corpus_path}/{OBS_ARRAY_NAME}")


def transform(obs: pd.DataFrame) -> pd.DataFrame:
    df = (
        obs.df[:]
        .groupby(
            by=[
                "dataset_id",
                "cell_type_ontology_term_id",
                "tissue_ontology_term_id",
                "assay_ontology_term_id",
                "development_stage_ontology_term_id",
                "disease_ontology_term_id",
                "ethnicity_ontology_term_id",
                "sex_ontology_term_id",
                "organism_ontology_term_id",
            ],
            as_index=False,
        )
        .size()
    ).rename(columns={"size": "n_cells"})
    return df


def load(corpus_path: str, df: pd.DataFrame):
    uri = f"{corpus_path}/{CELL_COUNTS_CUBE_NAME}"
    create_empty_cube(uri, cell_counts_schema)
    tiledb.from_pandas(uri, df, mode="append")
    return uri


def create_cell_count_cube(corpus_path: str):
    """
    Create cell count cube and write to disk
    """
    obs = extract(corpus_path)
    df = transform(obs)
    uri = load(corpus_path, df)
    logger.info(f"Cell count cube created and stored at {uri}")
