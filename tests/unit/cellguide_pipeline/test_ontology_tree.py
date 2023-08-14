import json
import unittest

from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


def compare_dicts(dict1, dict2):
    """
    This function recursive compares two dictionaries handling cases where the values
    are unordered arrays with elements that could be dictionaries.
    """
    if len(dict1) != len(dict2):
        return False

    for key in dict1:
        if key not in dict2:
            return False

        value1 = dict1[key]
        value2 = dict2[key]

        if isinstance(value1, dict) and isinstance(value2, dict):
            if not compare_dicts(value1, value2):
                return False
        elif isinstance(value1, list) and isinstance(value2, list):
            if len(value1) != len(value2):
                return False
            # check if the lists contain dictionaries as elements
            if len(value1) > 0 and isinstance(value1[0], dict) and isinstance(value2[0], dict):
                for i in range(len(value1)):
                    if not compare_dicts(value1[i], value2[i]):
                        return False
            elif sorted(value1) != sorted(value2):
                return False
        else:
            if value1 != value2:
                return False

    return True


class OntologyTreeBuilderTests(unittest.TestCase):
    def test__ontology_tree_builder(self):
        with open("tests/unit/cellguide_pipeline/fixtures/ontology_graph.json", "r") as f:
            expected__ontology_graph = json.load(f)
        with open("tests/unit/cellguide_pipeline/fixtures/all_states_per_cell_type.json", "r") as f:
            expected__all_states_per_cell_type = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            cell_counts_df = snapshot.cell_counts_cube.df[:]
            tree_builder = OntologyTreeBuilder(cell_counts_df)

            ontology_graph = tree_builder.ontology_graph
            self.assertTrue(compare_dicts(ontology_graph, expected__ontology_graph))

            all_states_per_cell_type = tree_builder.get_ontology_tree_state_per_celltype()
            self.assertTrue(compare_dicts(all_states_per_cell_type, expected__all_states_per_cell_type))
