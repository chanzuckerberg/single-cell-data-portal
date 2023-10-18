import json
import unittest

from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils.dict_compare import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)
from tests.unit.cellguide_pipeline.constants import (
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    CELLTYPE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME,
    ONTOLOGY_GRAPH_FIXTURE_FILENAME,
    TISSUE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class OntologyTreeBuilderTests(unittest.TestCase):
    def test__ontology_tree_builder(self):
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ONTOLOGY_GRAPH_FIXTURE_FILENAME}", "r") as f:
            expected__ontology_graph = json.load(f)
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{CELLTYPE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME}", "r") as f:
            expected__all_states_per_cell_type = json.load(f)
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{TISSUE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME}", "r") as f:
            expected__all_states_per_tissue = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            cell_counts_df = snapshot.cell_counts_cube.df[:]
            tree_builder = OntologyTreeBuilder(cell_counts_df)

            ontology_graph = convert_dataclass_to_dict_and_strip_nones(tree_builder.get_ontology_tree())
            self.assertTrue(compare_dicts(ontology_graph, expected__ontology_graph))

            all_states_per_cell_type = convert_dataclass_to_dict_and_strip_nones(
                tree_builder.get_ontology_tree_state_per_celltype()
            )
            self.assertTrue(compare_dicts(all_states_per_cell_type, expected__all_states_per_cell_type))

            all_states_per_tissue = convert_dataclass_to_dict_and_strip_nones(
                tree_builder.get_ontology_tree_state_per_tissue()
            )
            self.assertTrue(compare_dicts(all_states_per_tissue, expected__all_states_per_tissue))
