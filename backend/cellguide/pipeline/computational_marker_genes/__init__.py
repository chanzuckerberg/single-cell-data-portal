import logging

from backend.cellguide.pipeline.computational_marker_genes.computational_markers import MarkerGenesCalculator
from backend.cellguide.pipeline.constants import COMPUTATIONAL_MARKER_GENES_FOLDERNAME
from backend.cellguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import output_json_per_key
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot

logger = logging.getLogger(__name__)


def run(*, output_directory: str):
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    ontology_tree = get_ontology_tree_builder(snapshot=snapshot)
    marker_genes = get_computational_marker_genes(
        snapshot=snapshot,
        ontology_tree=ontology_tree,
    )
    output_json_per_key(marker_genes, f"{output_directory}/{COMPUTATIONAL_MARKER_GENES_FOLDERNAME}")


def get_computational_marker_genes(*, snapshot: WmgSnapshot, ontology_tree: OntologyTreeBuilder) -> dict:
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
    all_cell_types_in_corpus = ontology_tree.all_cell_type_ids_in_corpus

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
