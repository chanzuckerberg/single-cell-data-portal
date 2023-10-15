import json
import logging
import os

import tiledb

from backend.wmg.data.ontology_labels import gene_term_label, ontology_term_label
from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME,
    PRIMARY_FILTER_DIMENSIONS_FILENAME,
)
from backend.wmg.data.tissue_mapper import TissueMapper
from backend.wmg.data.utils import log_func_runtime
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
    PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG,
)
from backend.wmg.pipeline.utils import (
    load_pipeline_state,
    write_pipeline_state,
)

logger = logging.getLogger(__name__)


@log_func_runtime
def create_primary_filter_dimensions(*, corpus_path: str):
    """
    This method creates the primary filter dimensions for the WMG snapshot.
    """
    pipeline_state = load_pipeline_state(corpus_path)
    if not pipeline_state.get(EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG):
        raise ValueError(
            "'expression_summary_default' array does not exist. Please run 'create_expression_summary_default_cube' first."
        )
    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise ValueError(
            "'cell_counts' array does not exist. Please run 'create_cell_counts_cube_and_filter_relationships' first."
        )

    snapshot_id = os.path.basename(os.path.normpath(corpus_path))

    logger.info("Creating primary filter dimensions")
    with (
        tiledb.open(f"{corpus_path}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}", "r") as expression_summary_cube,
        tiledb.open(f"{corpus_path}/{CELL_COUNTS_CUBE_NAME}", "r") as cell_counts_cube,
    ):
        # get dataframes
        cell_counts_df = cell_counts_cube.df[:]
        expr_df = expression_summary_cube.df[:]

        # genes
        organism_gene_ids = list_grouped_primary_filter_dimensions_term_ids(
            expr_df, "gene_ontology_term_id", "organism_ontology_term_id"
        )
        organism_gene_terms = {
            organism_term_id: [{g: gene_term_label(g)} for g in gene_term_ids]
            for organism_term_id, gene_term_ids in organism_gene_ids.items()
        }

        # tissues
        organism_tissue_ids = list_grouped_primary_filter_dimensions_term_ids(
            cell_counts_df, "tissue_ontology_term_id", group_by_dim="organism_ontology_term_id"
        )

        # organisms
        organism_tissue_terms = {
            organism_term_id: [{t: ontology_term_label(t)} for t in order_tissues(tissue_term_ids)]
            for organism_term_id, tissue_term_ids in organism_tissue_ids.items()
        }

        # collate
        result = dict(
            snapshot_id=str(snapshot_id),
            organism_terms=[
                {o: ontology_term_label(o)} for o in sorted(cell_counts_df["organism_ontology_term_id"].unique())
            ],
            tissue_terms=organism_tissue_terms,
            gene_terms=organism_gene_terms,
        )
        logger.info("Writing primary filter dimensions")
        with open(f"{corpus_path}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}", "w") as f:
            json.dump(result, f)

        pipeline_state[PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG] = True
        write_pipeline_state(pipeline_state, corpus_path)


def list_grouped_primary_filter_dimensions_term_ids(
    df, primary_dim_name: str, group_by_dim: str
) -> dict[str, list[str]]:
    """
    This function takes a dataframe and two dimension names as input. It groups the dataframe by the second dimension name,
    and lists the unique values of the first dimension name for each group. The output is a dictionary where the keys are the
    unique values of the second dimension name, and the values are lists of unique values of the first dimension name for each group.

    :param df: The input dataframe.
    :param primary_dim_name: The name of the first dimension.
    :param group_by_dim: The name of the second dimension to group by.
    :return: A dictionary of lists of unique values of the first dimension for each group.
    """
    return (
        df[[primary_dim_name, group_by_dim]]
        .drop_duplicates()
        .groupby(group_by_dim)
        .agg(list)
        .to_dict()[primary_dim_name]
    )


def order_tissues(ontology_term_ids: list[str]) -> list[str]:
    """
    Order tissues based on appearance in TissueMapper.HIGH_LEVEL_TISSUES. This will maintain the priority set in
    that class which is intended to keep most relevant tissues on top and tissues that are related to be placed
    sequentially
    """
    ontology_term_ids = set(ontology_term_ids)
    ordered_ontology_term_ids = []
    for tissue in TissueMapper.HIGH_LEVEL_TISSUES:
        tissue = TissueMapper.reformat_ontology_term_id(tissue, to_writable=True)
        if tissue in ontology_term_ids:
            ontology_term_ids.remove(tissue)
            ordered_ontology_term_ids.append(tissue)

    if ontology_term_ids:
        ordered_ontology_term_ids = ordered_ontology_term_ids + list(ontology_term_ids)

    return ordered_ontology_term_ids
