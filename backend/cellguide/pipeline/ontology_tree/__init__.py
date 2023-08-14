from backend.cellguide.pipeline.constants import ONTOLOGY_TREE_FILENAME
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import output_json
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot


def run(output_directory):
    snapshot: WmgSnapshot = load_snapshot()
    cell_counts_df = snapshot.cell_counts_cube.df[:]

    tree_builder = OntologyTreeBuilder(cell_counts_df)
    output_json(tree_builder.ontology_graph, f"{output_directory}/{ONTOLOGY_TREE_FILENAME}")
