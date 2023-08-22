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


def run(*, output_directory) -> OntologyTreeBuilder:
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    ontology_graph, ontology_state_cell_types, ontology_state_tissues = get_ontology_tree_data(snapshot=snapshot)

    output_json(ontology_graph, f"{output_directory}/{ONTOLOGY_TREE_FILENAME}")
    output_json_per_key(ontology_state_cell_types, f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_CELLTYPE_FOLDERNAME}")
    output_json(ontology_state_tissues, f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_TISSUE_FOLDERNAME}")


def get_ontology_tree_builder(*, snapshot: WmgSnapshot):
    cell_counts_df = snapshot.cell_counts_cube.df[:]
    return OntologyTreeBuilder(cell_counts_df)


def get_ontology_tree_data(*, snapshot: WmgSnapshot) -> OntologyTreeData:
    tree_builder = get_ontology_tree_builder(snapshot=snapshot)
    ontology_graph = tree_builder.get_ontology_tree()
    all_states_per_cell_type = tree_builder.get_ontology_tree_state_per_celltype()
    all_states_per_tissue = tree_builder.get_ontology_tree_state_per_tissue()

    return OntologyTreeData(
        **{
            "tree_builder": tree_builder,
            "ontology_graph": ontology_graph,
            "all_states_per_cell_type": all_states_per_cell_type,
            "all_states_per_tissue": all_states_per_tissue,
        }
    )
