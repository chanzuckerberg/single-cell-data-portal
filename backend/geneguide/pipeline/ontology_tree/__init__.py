from backend.geneguide.pipeline.constants import (
    ONTOLOGY_TREE_FILENAME,
    ONTOLOGY_TREE_STATE_PER_GOTERM_FOLDERNAME,
    ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME,
)
from backend.geneguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.geneguide.pipeline.ontology_tree.types import OntologyTreeData
from backend.geneguide.pipeline.utils import output_json, output_json_per_key


def run(*, output_directory):
    ontology_tree_data = get_ontology_tree_data()

    output_json(
        ontology_tree_data.ontology_graph,
        f"{output_directory}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{ONTOLOGY_TREE_FILENAME}",
    )
    output_json_per_key(
        ontology_tree_data.all_states_per_go_term,
        f"{output_directory}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{ONTOLOGY_TREE_STATE_PER_GOTERM_FOLDERNAME}",
    )


def get_ontology_tree_builder() -> OntologyTreeBuilder:
    return OntologyTreeBuilder()


def get_ontology_tree_data() -> OntologyTreeData:
    tree_builder = get_ontology_tree_builder()
    ontology_graph = tree_builder.get_ontology_tree()
    all_states_per_go_term = tree_builder.get_ontology_tree_state_per_goterm()

    return OntologyTreeData(
        ontology_graph=ontology_graph,
        all_states_per_go_term=all_states_per_go_term,
    )
