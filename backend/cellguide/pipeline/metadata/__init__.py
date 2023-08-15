from backend.cellguide.pipeline.metadata.metadata_generator import generate_cellguide_card_metadata
from backend.cellguide.pipeline.metadata.types import CellMetadata
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder


def run(ontology_tree: OntologyTreeBuilder) -> dict[str, CellMetadata]:
    return generate_cellguide_card_metadata(ontology_tree.all_cell_type_ids_in_corpus)
