from backend.geneguide.pipeline.constants import GENEGUIDE_METADATA_FILENAME
from backend.geneguide.pipeline.metadata.metadata_generator import (
    generate_card_metadata,
)
from backend.geneguide.pipeline.metadata.types import GoTermMetadata
from backend.geneguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.geneguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.geneguide.pipeline.utils import output_json


def run(*, output_directory: str):
    ontology_tree = get_ontology_tree_builder()
    go_metadata = get_go_metadata(ontology_tree=ontology_tree)

    output_json(go_metadata, f"{output_directory}/{GENEGUIDE_METADATA_FILENAME}")


def get_go_metadata(*, ontology_tree: OntologyTreeBuilder) -> dict[str, GoTermMetadata]:
    return generate_card_metadata(ontology_tree.all_go_ids)
