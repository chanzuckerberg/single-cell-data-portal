import json
import os
import unittest

from backend.cellguide.pipeline.ontology_tree import get_celltype_to_tissue_mapping
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils import compare_dicts
from tests.unit.backend.cellguide.pipeline.constants import (
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    CELLTYPE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME,
    CELLTYPE_TO_TISSUE_MAPPING_FILENAME,
    ONTOLOGY_GRAPH_FIXTURE_FILENAME,
    ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME,
    TISSUE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class OntologyTreeBuilderTests(unittest.TestCase):
    def test__ontology_tree_builder(self):
        organisms = next(os.walk(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/"))[1]

        for organism in organisms:
            with open(
                f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{organism}/{ONTOLOGY_GRAPH_FIXTURE_FILENAME}",
                "r",
            ) as f:
                expected__ontology_graph = json.load(f)
            with open(
                f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{organism}/{CELLTYPE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME}",
                "r",
            ) as f:
                expected__all_states_per_cell_type = json.load(f)
            with open(
                f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{organism}/{TISSUE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME}",
                "r",
            ) as f:
                expected__all_states_per_tissue = json.load(f)
            with open(
                f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{organism}/{CELLTYPE_TO_TISSUE_MAPPING_FILENAME}",
                "r",
            ) as f:
                expected__celltype_to_tissue_mapping = json.load(f)
            with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
                cell_counts_df = snapshot.cell_counts_cube.df[:]
                cell_counts_df = cell_counts_df[
                    cell_counts_df["organism_ontology_term_id"] == organism.replace("_", ":")
                ]
                tree_builder = OntologyTreeBuilder(cell_counts_df)

                ontology_graph = convert_dataclass_to_dict_and_strip_nones(tree_builder.get_ontology_tree())
                self.assertTrue(compare_dicts(ontology_graph, expected__ontology_graph))

                all_states_per_cell_type = convert_dataclass_to_dict_and_strip_nones(
                    tree_builder.get_ontology_tree_state_per_celltype()
                )
                self.assertTrue(compare_dicts(all_states_per_cell_type, expected__all_states_per_cell_type))

                all_states_per_tissue = tree_builder.get_ontology_tree_state_per_tissue()
                self.assertTrue(
                    compare_dicts(
                        convert_dataclass_to_dict_and_strip_nones(all_states_per_tissue),
                        expected__all_states_per_tissue,
                    )
                )

                celltype_to_tissue_mapping = get_celltype_to_tissue_mapping(all_states_per_tissue)
                self.assertTrue(compare_dicts(celltype_to_tissue_mapping, expected__celltype_to_tissue_mapping))
