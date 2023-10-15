import json
import logging
import os

import pandas as pd
import tiledb

from backend.cellguide.pipeline.computational_marker_genes.computational_markers import MarkerGenesCalculator
from backend.cellguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.wmg.data.schemas.marker_gene_cube_schema import marker_genes_schema
from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME,
    MARKER_GENES_CUBE_NAME,
    PRIMARY_FILTER_DIMENSIONS_FILENAME,
    WmgSnapshot,
)
from backend.wmg.data.utils import log_func_runtime
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
    MARKER_GENES_CUBE_CREATED_FLAG,
    PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG,
)
from backend.wmg.pipeline.utils import (
    create_empty_cube_if_needed,
    load_pipeline_state,
    write_pipeline_state,
)

logger = logging.getLogger(__name__)


@log_func_runtime
def create_marker_genes_cube(*, corpus_path: str):
    pipeline_state = load_pipeline_state(corpus_path)

    expression_summary_default_cube_uri = os.path.join(corpus_path, EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME)
    cell_counts_cube_uri = os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)
    if not pipeline_state.get(EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG):
        raise ValueError(
            "'expression_summary_default' cube does not exist. Please run 'create_expression_summary_default_cube' first."
        )

    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise ValueError(
            "'cell_counts' cube does not exist. Please run 'create_cell_counts_cube_and_filter_relationships' first."
        )

    if not pipeline_state.get(PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG):
        raise ValueError(
            "'primary_filter_dimensions' file does not exist. Please run 'create_primary_filter_dimensions' first."
        )

    logger.info("Calculating marker genes.")
    with (
        open(os.path.join(corpus_path, PRIMARY_FILTER_DIMENSIONS_FILENAME), "r") as f,
        tiledb.open(cell_counts_cube_uri, "r") as cell_counts_cube,
        tiledb.open(expression_summary_default_cube_uri, "r") as expression_summary_default_cube,
    ):
        primary_filter_dimensions = json.load(f)
        snapshot = WmgSnapshot(
            primary_filter_dimensions=primary_filter_dimensions,
            cell_counts_cube=cell_counts_cube,
            expression_summary_default_cube=expression_summary_default_cube,
        )
        ontology_tree = get_ontology_tree_builder(snapshot=snapshot)
        calculator = MarkerGenesCalculator(
            snapshot=snapshot,
            all_cell_type_ids_in_corpus=ontology_tree.all_cell_type_ids_in_corpus,
            groupby_terms=["organism_ontology_term_id", "tissue_ontology_term_id"],
        )
        marker_genes = calculator.get_computational_marker_genes()
    marker_gene_records = []
    for key in marker_genes:
        for marker_gene in marker_genes[key]:
            marker_gene_records.append(
                {
                    "tissue_ontology_term_id": marker_gene.groupby_dims["tissue_ontology_term_id"],
                    "organism_ontology_term_id": marker_gene.groupby_dims["organism_ontology_term_id"],
                    "cell_type_ontology_term_id": key,
                    "gene_ontology_term_id": marker_gene.gene_ontology_term_id,
                    "marker_score": marker_gene.marker_score,
                    "specificity": marker_gene.specificity,
                }
            )
    marker_genes_df = pd.DataFrame(marker_gene_records)
    uri = os.path.join(corpus_path, MARKER_GENES_CUBE_NAME)
    create_empty_cube_if_needed(uri, marker_genes_schema)
    logger.info("Writing marker genes cube.")
    tiledb.from_pandas(uri, marker_genes_df, mode="append")

    pipeline_state[MARKER_GENES_CUBE_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
