import logging

from backend.cellguide.pipeline.computational_marker_genes.computational_markers import MarkerGenesCalculator
from backend.cellguide.pipeline.constants import COMPUTATIONAL_MARKER_GENES_FOLDERNAME, MARKER_GENE_PRESENCE_FILENAME
from backend.cellguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import output_json, output_json_per_key
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot

logger = logging.getLogger(__name__)


def run(*, output_directory: str):
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    ontology_tree = get_ontology_tree_builder(snapshot=snapshot)
    marker_genes, reformatted_marker_genes = get_computational_marker_genes(
        snapshot=snapshot,
        ontology_tree=ontology_tree,
    )
    output_json_per_key(marker_genes, f"{output_directory}/{COMPUTATIONAL_MARKER_GENES_FOLDERNAME}")
    output_json(
        marker_genes, f"{output_directory}/{COMPUTATIONAL_MARKER_GENES_FOLDERNAME}/{MARKER_GENE_PRESENCE_FILENAME}"
    )


def get_computational_marker_genes(*, snapshot: WmgSnapshot, ontology_tree: OntologyTreeBuilder) -> tuple[dict, dict]:
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

    # gene --> organism --> tissue --> cell types
    reformatted_marker_genes = {}
    for cell_type_id, marker_gene_stats_list in marker_genes.items():
        for marker_gene_stats in marker_gene_stats_list:
            symbol = marker_gene_stats["symbol"]
            tissue = marker_gene_stats["groupby_dims"].get("tissue_ontology_term_label", "All Tissues")
            organism = marker_gene_stats["groupby_dims"]["organism_ontology_term_label"]

            if symbol not in reformatted_marker_genes:
                reformatted_marker_genes[symbol] = {}
            if organism not in reformatted_marker_genes[symbol]:
                reformatted_marker_genes[symbol][organism] = {}
            if tissue not in reformatted_marker_genes[symbol][organism]:
                reformatted_marker_genes[symbol][organism][tissue] = []
            reformatted_marker_genes[symbol][organism][tissue].append(cell_type_id)

    # make the lists unique in reformatted_marker_genes
    for symbol in reformatted_marker_genes:
        for organism in reformatted_marker_genes[symbol]:
            for tissue in reformatted_marker_genes[symbol][organism]:
                reformatted_marker_genes[symbol][organism][tissue] = list(
                    set(reformatted_marker_genes[symbol][organism][tissue])
                )

    return marker_genes, reformatted_marker_genes
