import logging
import time

import tiledb

from backend.wmg.data.schemas.cube_schema import cell_counts_schema, expression_summary_schema
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME, EXPRESSION_SUMMARY_CUBE_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.wmg_cube import logger
from backend.atlas_asset_pipelines.cubes.load import load_data_into_cube

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def create_empty_cube(uri: str, schema):
    """
    Create an empty cube with expected schema (dimensions and attributes) at given uri
    """
    tiledb.Array.create(uri, schema, overwrite=True)


def create_cubes(corpus_path):
    create_expression_summary_cube(corpus_path)
    create_cell_count_cube(corpus_path)


def create_cell_count_cube(corpus_path: str):
    """
    Create cell count cube and write to disk
    """
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


def create_expression_summary_cube(corpus_path):
    """
    Create queryable cube and write to disk
    """
    uri = f"{corpus_path}/{EXPRESSION_SUMMARY_CUBE_NAME}"
    cube_creation_start_time = time.time()

    with tiledb.scope_ctx(create_ctx()):
        # Create cube
        create_empty_cube(uri, expression_summary_schema)

        # load data
        dims, vals = load_data_into_cube(corpus_path, uri)

        logger.debug("Saving cube to tiledb")
        with tiledb.open(uri, "w") as cube:
            cube[tuple(dims)] = vals

        logger.debug("Cube created, start consolidation")
        tiledb.consolidate(uri)

        logger.debug("Cube consolidated, start vacuumming")
        tiledb.vacuum(uri)

    logger.debug("Cube creation complete")
    create_cube_sec = time.time() - cube_creation_start_time
    logger.info(f"Expression summary cube: time to create {create_cube_sec}, uri={uri}")



