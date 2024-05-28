import json
import unittest

from backend.cellguide.pipeline.metadata.metadata_generator import (
    generate_cellguide_card_metadata,
    generate_cellguide_tissue_card_metadata,
)
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from backend.common.census_cube.utils import get_all_cell_type_ids_in_corpus, get_all_tissue_ids_in_corpus
from tests.test_utils import compare_dicts
from tests.unit.backend.cellguide.pipeline.constants import (
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    CELLTYPE_METADATA_FIXTURE_FILENAME,
    TISSUE_METADATA_FIXTURE_FILENAME,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class TestMetadataGenerator(unittest.TestCase):
    def test__cell_metadata_generator(self):
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{CELLTYPE_METADATA_FIXTURE_FILENAME}", "r") as f:
            expected__cell_metadata = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            all_cell_type_ids_in_corpus = get_all_cell_type_ids_in_corpus(snapshot=snapshot)
            cell_metadata = generate_cellguide_card_metadata(all_cell_type_ids_in_corpus)
            self.assertTrue(
                compare_dicts(convert_dataclass_to_dict_and_strip_nones(cell_metadata), expected__cell_metadata)
            )

    def test__tissue_metadata_generator(self):
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{TISSUE_METADATA_FIXTURE_FILENAME}", "r") as f:
            expected__tissue_metadata = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            all_tissue_ids_in_corpus = get_all_tissue_ids_in_corpus(snapshot=snapshot)
            tissue_metadata = generate_cellguide_tissue_card_metadata(all_tissue_ids_in_corpus=all_tissue_ids_in_corpus)
            self.assertTrue(
                compare_dicts(convert_dataclass_to_dict_and_strip_nones(tissue_metadata), expected__tissue_metadata)
            )
