from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME

import logging

import tiledb

from backend.wmg.data.schemas.cube_schema import cell_counts_schema
from backend.wmg.data.utils import create_empty_cube

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def create_cell_count_cube(corpus_path: str):
    """
    Create cell count cube and write to disk
    """
    logger.info("Creating cell count cube")
    uri = f"{corpus_path}/{CELL_COUNTS_CUBE_NAME}"
    with tiledb.open(f"{corpus_path}/obs") as obs:
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
        )
        df = df.rename(columns={"size": "n_cells"})
        create_empty_cube(uri, cell_counts_schema)
        tiledb.from_pandas(uri, df, mode="append")
        logger.info("Cell count cube creation complete")
