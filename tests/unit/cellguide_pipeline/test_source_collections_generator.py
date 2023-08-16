import json
import unittest

from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.source_collections.source_collections_generator import generate_source_collections_data
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils.dict_compare import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class TestSourceCollectionsGenerator(unittest.TestCase):
    def test__source_collections_generator(self):
        with open("tests/unit/cellguide_pipeline/fixtures/source_collections.json", "r") as f:
            expected__source_collections = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            cell_counts_df = snapshot.cell_counts_cube.df[:]
            tree_builder = OntologyTreeBuilder(cell_counts_df)

            source_collections = generate_source_collections_data(
                all_cell_type_ids_in_corpus=tree_builder.all_cell_type_ids_in_corpus
            )
            self.assertTrue(
                compare_dicts(
                    convert_dataclass_to_dict_and_strip_nones(source_collections), expected__source_collections
                )
            )
