from backend.corpus_asset_pipelines.expression_summary_cube.extract import extract_var_data
from backend.corpus_asset_pipelines.expression_summary_cube.job import create_expression_summary_cube
from backend.corpus_asset_pipelines.expression_summary_cube.load import build_in_mem_cube
from backend.corpus_asset_pipelines.expression_summary_cube.transform import cube_indexed_dims_no_gene_ontology, \
    reduce_X, make_cube_index
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME

import logging
import time

import numpy as np
import tiledb

from backend.wmg.data.schemas.cube_schema import cube_non_indexed_dims, cell_counts_schema
from backend.wmg.data.tiledb import create_ctx
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


def create_cubes(corpus_path):
    create_expression_summary_cube(corpus_path)
    create_cell_count_cube(corpus_path)


def load_data_into_cube(tdb_group, uri: str):
    """
    Load data from the concatenated corpus into the queryable cube
    """
    ctx = create_ctx()
    start_time = time.time()
    logger.debug(f"Start loading big cube at : {uri}")

    gene_ontology_term_ids = extract_var_data(tdb_group, ctx)
    n_genes = len(gene_ontology_term_ids)

    ##
    # Reduce X
    ##
    big_cube_atts = cube_indexed_dims_no_gene_ontology + cube_non_indexed_dims
    cell_labels, cube_index = make_cube_index(tdb_group, big_cube_atts)
    n_groups = len(cube_index)

    cube_sum = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_nnz = np.zeros((n_groups, n_genes), dtype=np.uint64)
    cube_min = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_max = np.zeros((n_groups, n_genes), dtype=np.float32)

    # pass 1 - sum, nnz, min, max
    reduce_X(tdb_group, start_time, cell_labels.cube_idx.values, cube_sum, cube_nnz, cube_min, cube_max)

    return build_in_mem_cube(gene_ontology_term_ids, cube_index, cube_non_indexed_dims, cube_sum, cube_nnz)




