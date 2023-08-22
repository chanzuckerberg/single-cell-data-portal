import logging

from backend.cellguide.pipeline.computational_marker_genes.computational_markers import MarkerGenesCalculator
from backend.cellguide.pipeline.computational_marker_genes.types import ComputationalMarkerGenes
from backend.cellguide.pipeline.constants import COMPUTATIONAL_MARKER_GENES_FOLDERNAME
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import output_json
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot

logger = logging.getLogger(__name__)


def run(*, output_directory: str, ontology_tree: OntologyTreeBuilder):
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    marker_genes = get_marker_genes_per_and_across_tissues(
        snapshot=snapshot,
        all_cell_types_in_corpus=ontology_tree.all_cell_type_ids_in_corpus,
    )
    output_marker_genes(marker_genes, output_directory)


def get_marker_genes_per_and_across_tissues(*, snapshot: WmgSnapshot, all_cell_types_in_corpus: list[str]) -> dict:
    """
    This function calculates the marker genes per tissue and across tissues.

    Arguments
    ---------
    snapshot - The WMG snapshot
    all_cell_types_in_corpus - List of all cell types in the corpus.

    Returns
    -------
    dict - A dictionary containing the marker genes per tissue and across tissues keyed by cell type ontology term ID.
    """

    calculator = MarkerGenesCalculator(
        snapshot=snapshot,
        all_cell_type_ids_in_corpus=all_cell_types_in_corpus,
        groupby_terms=["organism_ontology_term_id", "cell_type_ontology_term_id"],
    )
    marker_genes = calculator.get_computational_marker_genes()

    calculator = MarkerGenesCalculator(
        snapshot=snapshot,
        all_cell_type_ids_in_corpus=all_cell_types_in_corpus,
        groupby_terms=["organism_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"],
    )
    marker_genes_per_tissue = calculator.get_computational_marker_genes()

    for key in marker_genes_per_tissue:
        if key in marker_genes:
            marker_genes[key] += marker_genes_per_tissue[key]
        else:
            marker_genes[key] = marker_genes_per_tissue[key]

    return marker_genes


def output_marker_genes(marker_genes: dict[str, list[ComputationalMarkerGenes]], output_directory: str):
    """
    This function outputs the marker genes to a specified directory. Each cell type will be output to its own JSON file.

    Arguments
    ---------
    marker_genes - A dictionary containing the marker genes per tissue and across tissues keyed by cell type ontology term ID.
    output_directory - The directory where the output files will be saved.
    """

    for key in marker_genes:
        output_json(marker_genes[key], f"{output_directory}/{COMPUTATIONAL_MARKER_GENES_FOLDERNAME}/{key}.json")
