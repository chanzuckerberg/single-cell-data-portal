from backend.cellguide.pipeline.constants import (
    ONTOLOGY_TREE_FILENAME,
    ONTOLOGY_TREE_STATE_PER_CELLTYPE_FOLDERNAME,
    ONTOLOGY_TREE_STATE_PER_TISSUE_FOLDERNAME,
)
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.ontology_tree.types import OntologyTreeData
from backend.cellguide.pipeline.utils import output_json, output_json_per_key
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot


def run(*, output_directory):
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    _, ontology_tree_data = get_ontology_tree_data(snapshot=snapshot)

    output_json(ontology_tree_data.ontology_graph, f"{output_directory}/{ONTOLOGY_TREE_FILENAME}")
    output_json_per_key(
        ontology_tree_data.all_states_per_cell_type, f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_CELLTYPE_FOLDERNAME}"
    )
    output_json_per_key(
        ontology_tree_data.all_states_per_tissue, f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_TISSUE_FOLDERNAME}"
    )


def get_ontology_tree_builder(*, snapshot: WmgSnapshot) -> OntologyTreeBuilder:
    cell_counts_df = snapshot.cell_counts_cube.df[:]
    return OntologyTreeBuilder(cell_counts_df)


def get_ontology_tree_data(*, snapshot: WmgSnapshot) -> tuple[OntologyTreeBuilder, OntologyTreeData]:
    tree_builder = get_ontology_tree_builder(snapshot=snapshot)
    ontology_graph = tree_builder.get_ontology_tree()
    all_states_per_cell_type = tree_builder.get_ontology_tree_state_per_celltype()
    all_states_per_tissue = tree_builder.get_ontology_tree_state_per_tissue()

    return tree_builder, OntologyTreeData(
        ontology_graph=ontology_graph,
        all_states_per_cell_type=all_states_per_cell_type,
        all_states_per_tissue=all_states_per_tissue,
    )
