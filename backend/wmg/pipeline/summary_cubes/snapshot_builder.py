import logging
import os
import pathlib

import cellxgene_census
import tiledbsoma as soma

from backend.wmg.data.constants import (
    GENE_EXPRESSION_COUNT_MIN_THRESHOLD,
)
from backend.wmg.data.utils import log_func_runtime
from backend.wmg.pipeline.summary_cubes.cell_counts_and_filters import create_cell_counts_cube_and_filter_relationships
from backend.wmg.pipeline.summary_cubes.dataset_metadata import create_dataset_metadata
from backend.wmg.pipeline.summary_cubes.expression_summary import ExpressionSummaryCubeBuilder
from backend.wmg.pipeline.summary_cubes.expression_summary_default import create_expression_summary_default_cube
from backend.wmg.pipeline.summary_cubes.marker_genes import create_marker_genes_cube
from backend.wmg.pipeline.summary_cubes.primary_filter_dimensions import create_primary_filter_dimensions

logger = logging.getLogger(__name__)


class WmgSnapshotBuilder:
    def __init__(self, *, corpus_path: str, organismInfo: dict):
        """
        Args:
            corpus_path (str): The path to the corpus.
            organismInfo (dict): Information about the organism.
                ex:
                organismInfo = {
                    "id": "NCBITaxon:9606",
                    "label": "homo_sapiens",
                }
                Note that "label" should be an organism label used by census.
        """

        self.snapshot_id = os.path.basename(os.path.normpath(corpus_path))
        self.organism = organismInfo["label"]
        self.organismId = organismInfo["id"]
        self.corpus_path = os.path.join(corpus_path, organismInfo["id"].replace(":", "_"))

        logger.info(f"Creating directory {self.corpus_path}")
        pathlib.Path(self.corpus_path).mkdir(parents=True, exist_ok=True)

    @log_func_runtime
    def run_pipeline(self):
        """
        This method runs the entire pipeline for creating summary cubes. It sequentially calls the methods for creating
        the expression summary cube, the default expression summary cube, the cell counts cube and filter relationships,
        the primary filter dimensions, and the marker genes cube.
        """
        with cellxgene_census.open_soma(
            census_version="latest",
        ) as census:
            organism = census["census_data"][self.organism]
            with organism.axis_query(
                "RNA",
                obs_query=soma.AxisQuery(
                    value_filter=f"is_primary_data == True and nnz >= {GENE_EXPRESSION_COUNT_MIN_THRESHOLD}"
                ),
            ) as query:
                ExpressionSummaryCubeBuilder(query=query, corpus_path=self.corpus_path).create_expression_summary_cube()
                create_expression_summary_default_cube(corpus_path=self.corpus_path)
                create_cell_counts_cube_and_filter_relationships(query=query, corpus_path=self.corpus_path)
                create_primary_filter_dimensions(
                    corpus_path=self.corpus_path, organismId=self.organismId, snapshot_id=self.snapshot_id
                )
                create_marker_genes_cube(corpus_path=self.corpus_path)
                create_dataset_metadata(corpus_path=self.corpus_path)
