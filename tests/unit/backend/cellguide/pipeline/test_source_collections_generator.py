import json
import unittest
from unittest.mock import patch

from backend.cellguide.pipeline.source_collections.source_collections_generator import generate_source_collections_data
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from backend.common.census_cube.utils import get_all_cell_type_ids_in_corpus
from tests.test_utils import compare_dicts
from tests.test_utils.mocks import (
    mock_get_collections_from_curation_endpoint,
    mock_get_datasets_from_curation_endpoint,
)
from tests.unit.backend.cellguide.pipeline.constants import (
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    SOURCE_COLLECTIONS_FIXTURE_FILENAME,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class TestSourceCollectionsGenerator(unittest.TestCase):
    def test__source_collections_generator(self):
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{SOURCE_COLLECTIONS_FIXTURE_FILENAME}", "r") as f:
            expected__source_collections = json.load(f)
        with (
            load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot,
            patch(
                "backend.cellguide.pipeline.source_collections.source_collections_generator.get_datasets_from_discover_api",
                new=mock_get_datasets_from_curation_endpoint,
            ),
            patch(
                "backend.cellguide.pipeline.source_collections.source_collections_generator.get_collections_from_discover_api",
                new=mock_get_collections_from_curation_endpoint,
            ),
        ):
            all_cell_type_ids_in_corpus = get_all_cell_type_ids_in_corpus(snapshot)
            source_collections = generate_source_collections_data(all_cell_type_ids_in_corpus)
            self.assertTrue(
                compare_dicts(
                    convert_dataclass_to_dict_and_strip_nones(source_collections), expected__source_collections
                )
            )
