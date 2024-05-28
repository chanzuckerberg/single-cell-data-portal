import logging

from backend.cellguide.common.constants import (
    COMPUTATIONAL_MARKER_GENES_FOLDERNAME,
    MARKER_GENE_DATA_FILENAME,
    MARKER_GENE_PRESENCE_FILENAME,
)
from backend.cellguide.common.data import format_marker_gene_data
from backend.cellguide.pipeline.constants import CELLGUIDE_CENSUS_CUBE_DATA_SCHEMA_VERSION
from backend.cellguide.pipeline.utils import output_json, output_json_per_key
from backend.common.census_cube.data import snapshot as sn
from backend.common.marker_genes.computational_markers import MarkerGenesCalculator
from backend.common.marker_genes.constants import MARKER_SCORE_THRESHOLD

logger = logging.getLogger(__name__)


def run(*, output_directory: str):
    snapshot = sn.load_snapshot(snapshot_schema_version=CELLGUIDE_CENSUS_CUBE_DATA_SCHEMA_VERSION)

    marker_genes, reformatted_marker_genes, formatted_marker_gene_data = get_computational_marker_genes(
        snapshot=snapshot,
    )
    output_json_per_key(marker_genes, f"{output_directory}/{COMPUTATIONAL_MARKER_GENES_FOLDERNAME}")
    output_json(
        reformatted_marker_genes,
        f"{output_directory}/{COMPUTATIONAL_MARKER_GENES_FOLDERNAME}/{MARKER_GENE_PRESENCE_FILENAME}",
    )
    output_json(
        formatted_marker_gene_data,
        f"{output_directory}/{COMPUTATIONAL_MARKER_GENES_FOLDERNAME}/{MARKER_GENE_DATA_FILENAME}",
    )


def get_computational_marker_genes(*, snapshot: sn.CensusCubeSnapshot) -> tuple[dict, dict]:
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
        groupby_terms=["organism_ontology_term_id", "cell_type_ontology_term_id"],
    )
    marker_genes = calculator.get_computational_marker_genes()

    calculator = MarkerGenesCalculator(
        snapshot=snapshot,
        groupby_terms=["organism_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"],
    )
    marker_genes_per_tissue = calculator.get_computational_marker_genes()

    for key in marker_genes_per_tissue:
        if key in marker_genes:
            marker_genes[key] += marker_genes_per_tissue[key]
        else:
            marker_genes[key] = marker_genes_per_tissue[key]

    # convert all groupby_dims IDs to labels as required by CellGuide
    organism_id_to_name = {k: v for d in snapshot.primary_filter_dimensions["organism_terms"] for k, v in d.items()}
    tissue_id_to_name = {
        k: v
        for organism in snapshot.primary_filter_dimensions["tissue_terms"]
        for i in snapshot.primary_filter_dimensions["tissue_terms"][organism]
        for k, v in i.items()
    }
    for _, marker_gene_stats_list in marker_genes.items():
        for marker_gene_stats in marker_gene_stats_list:
            groupby_dims = marker_gene_stats.groupby_dims
            groupby_terms = list(groupby_dims.keys())
            groupby_term_labels = [term.rsplit("_", 1)[0] + "_label" for term in groupby_terms]
            groupby_dims_new = dict(
                zip(groupby_term_labels, (groupby_dims[term] for term in groupby_terms), strict=False)
            )

            for key in groupby_dims_new:
                if key == "tissue_ontology_term_label":
                    groupby_dims_new[key] = tissue_id_to_name.get(groupby_dims_new[key], groupby_dims_new[key])
                elif key == "organism_ontology_term_label":
                    groupby_dims_new[key] = organism_id_to_name.get(groupby_dims_new[key], groupby_dims_new[key])

            marker_gene_stats.groupby_dims = groupby_dims_new

    reformatted_marker_genes = {}
    for cell_type_id, marker_gene_stats_list in marker_genes.items():
        for marker_gene_stats in marker_gene_stats_list:
            marker_score = marker_gene_stats.marker_score
            symbol = marker_gene_stats.symbol
            tissue = marker_gene_stats.groupby_dims.get("tissue_ontology_term_label", "All Tissues")
            organism = marker_gene_stats.groupby_dims["organism_ontology_term_label"]

            if marker_score <= MARKER_SCORE_THRESHOLD:
                continue

            if symbol not in reformatted_marker_genes:
                reformatted_marker_genes[symbol] = {}
            if organism not in reformatted_marker_genes[symbol]:
                reformatted_marker_genes[symbol][organism] = {}
            if tissue not in reformatted_marker_genes[symbol][organism]:
                reformatted_marker_genes[symbol][organism][tissue] = []

            data = dict(
                marker_score=marker_score,
                me=marker_gene_stats.me,
                pc=marker_gene_stats.pc,
                cell_type_id=cell_type_id,
            )
            reformatted_marker_genes[symbol][organism][tissue].append(data)

    # reformat the data to be a nested dictionary with structure organism-->tissue-->celltype-->genes
    organism_tissue_celltype_genes_data = format_marker_gene_data(reformatted_marker_genes)
    return marker_genes, reformatted_marker_genes, organism_tissue_celltype_genes_data
