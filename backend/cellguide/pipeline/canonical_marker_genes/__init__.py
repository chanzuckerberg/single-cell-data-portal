import logging

from backend.cellguide.pipeline.canonical_marker_genes.canonical_markers import CanonicalMarkerGenesCompiler
from backend.cellguide.pipeline.constants import CANONICAL_MARKER_GENES_FILENAME
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import output_json
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import load_snapshot

logger = logging.getLogger(__name__)


def run(*, output_directory: str, ontology_tree: OntologyTreeBuilder):
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    marker_gene_compiler = CanonicalMarkerGenesCompiler(snapshot.cell_counts_cube.df[:])
    parsed_asctb_table_entries = marker_gene_compiler.get_processed_asctb_table_entries()

    num_cell_types_in_corpus = len(
        set(parsed_asctb_table_entries).intersection(ontology_tree.all_cell_type_ids_in_corpus)
    )
    logger.info(
        f"Parsed {len(parsed_asctb_table_entries)} cell types in ASCTB, out of which {num_cell_types_in_corpus} are in CellGuide."
    )

    output_json(parsed_asctb_table_entries, f"{output_directory}/{CANONICAL_MARKER_GENES_FILENAME}")
