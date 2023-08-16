from backend.cellguide.pipeline.constants import (
    ONTOLOGY_TREE_FILENAME,
    ONTOLOGY_TREE_STATE_PER_CELLTYPE_FILENAME,
    ONTOLOGY_TREE_STATE_PER_TISSUE_FILENAME,
)
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import output_json
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot


def run(output_directory):
    snapshot: WmgSnapshot = load_snapshot()
    cell_counts_df = snapshot.cell_counts_cube.df[:]

    tree_builder = OntologyTreeBuilder(cell_counts_df)
    ontology_graph = tree_builder.get_ontology_tree()

    output_json(ontology_graph, f"{output_directory}/{ONTOLOGY_TREE_FILENAME}")

    all_states_per_cell_type = tree_builder.get_ontology_tree_state_per_celltype()
    output_json(all_states_per_cell_type, f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_CELLTYPE_FILENAME}")

    all_states_per_tissue = tree_builder.get_ontology_tree_state_per_tissue()
    output_json(all_states_per_tissue, f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_TISSUE_FILENAME}")
