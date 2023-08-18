from backend.cellguide.pipeline.constants import SOURCE_COLLECTIONS_FILENAME
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.source_collections.source_collections_generator import generate_source_collections_data
from backend.cellguide.pipeline.source_collections.types import SourceCollectionsData
from backend.cellguide.pipeline.utils import output_json


def run(*, output_directory: str, ontology_tree: OntologyTreeBuilder) -> dict[str, SourceCollectionsData]:
    data = generate_source_collections_data(ontology_tree.all_cell_type_ids_in_corpus)
    output_json(data, f"{output_directory}/{SOURCE_COLLECTIONS_FILENAME}")
    return data
