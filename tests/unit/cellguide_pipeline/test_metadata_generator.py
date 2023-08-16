import json
import unittest

from backend.cellguide.pipeline.metadata import run as run_metadata_pipeline
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.common.utils.dataclass import convert_dataclass_to_dict
from tests.test_utils.dict_compare import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class TestMetadataGenerator(unittest.TestCase):
    def test__cell_metadata_generator(self):
        with open("tests/unit/cellguide_pipeline/fixtures/cell_metadata.json", "r") as f:
            expected__cell_metadata = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            cell_counts_df = snapshot.cell_counts_cube.df[:]
            tree_builder = OntologyTreeBuilder(cell_counts_df)

            cell_metadata = run_metadata_pipeline(tree_builder)
            self.assertTrue(compare_dicts(convert_dataclass_to_dict(cell_metadata), expected__cell_metadata))
