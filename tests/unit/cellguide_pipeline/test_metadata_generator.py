import json
import unittest

from backend.cellguide.pipeline.metadata.metadata_generator import (
    generate_cellguide_card_metadata,
    generate_cellguide_tissue_card_metadata,
)
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
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

            cell_metadata = generate_cellguide_card_metadata(
                all_cell_type_ids_in_corpus=tree_builder.all_cell_type_ids_in_corpus
            )

            self.assertTrue(
                compare_dicts(convert_dataclass_to_dict_and_strip_nones(cell_metadata), expected__cell_metadata)
            )

    def test__tissue_metadata_generator(self):
        with open("tests/unit/cellguide_pipeline/fixtures/tissue_metadata.json", "r") as f:
            expected__tissue_metadata = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            cell_counts_df = snapshot.cell_counts_cube.df[:]
            tree_builder = OntologyTreeBuilder(cell_counts_df)

            tissue_metadata = generate_cellguide_tissue_card_metadata(
                all_tissue_ids_in_corpus=tree_builder.all_tissue_ids_in_corpus
            )
            self.assertTrue(
                compare_dicts(convert_dataclass_to_dict_and_strip_nones(tissue_metadata), expected__tissue_metadata)
            )
