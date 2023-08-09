from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot


def run():
    snapshot: WmgSnapshot = load_snapshot()
    cell_counts_df = snapshot.cell_counts_cube.df[:]

    tree_builder = OntologyTreeBuilder(cell_counts_df)
    tree_builder.get_ontology_tree_state_per_celltype()
    tree_builder.get_ontology_tree_state_per_tissue()

    # upload to S3 logic
    """... logic ..."""
