import os
import unittest
from unittest.mock import patch

from backend.cellguide.pipeline.gpt_descriptions.constants import (
    GPT_CELLTYPE_DESCRIPTION_SYSTEM_ROLE,
    GPT_CELLTYPE_DESCRIPTION_USER_ROLE,
    GPT_CELLTYPE_SEO_DESCRIPTION_USER_ROLE,
)
from backend.cellguide.pipeline.gpt_descriptions.gpt_description_generator import (
    generate_new_gpt_descriptions,
    generate_new_seo_gpt_descriptions,
)
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from tests.test_utils import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"

WHICH_OBJECTS_DO_NOT_EXIST = [
    {"id": "CL:0000576", "name": "monocyte", "fake_gpt_description": "fake description monocyte"},
    {"id": "CL:0005025", "name": "visceromotor neuron", "fake_gpt_description": "fake_description visceromotor neuron"},
]


class MockS3Provider:
    def does_object_exist(self, bucket_name: str, object_key: str):
        cell_type = os.path.basename(object_key).split(".json")[0]
        cell_type = cell_type.replace("_", ":")
        return cell_type not in [i["id"] for i in WHICH_OBJECTS_DO_NOT_EXIST]


class MockOpenAIProvider:
    def generate_gpt_output(self, *, system_role: str, user_role: str):
        return f"MOCK OUTPUT: {system_role} {user_role}"


class MockCellGuideConfig:
    def __init__(self):
        self.bucket = "cellguide-data-public-test"


class TestGptDescriptionGenerator(unittest.TestCase):
    def test__gpt_description_generator(self):
        expected_outputs = {}
        for obj in WHICH_OBJECTS_DO_NOT_EXIST:
            expected_outputs[
                obj["id"]
            ] = f"MOCK OUTPUT: {GPT_CELLTYPE_DESCRIPTION_SYSTEM_ROLE} {GPT_CELLTYPE_DESCRIPTION_USER_ROLE(obj['name'])}"

        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot, patch(
            "backend.cellguide.pipeline.gpt_descriptions.gpt_description_generator.S3Provider", new=MockS3Provider
        ), patch(
            "backend.cellguide.pipeline.gpt_descriptions.gpt_description_generator.OpenAIProvider",
            new=MockOpenAIProvider,
        ), patch(
            "backend.cellguide.pipeline.gpt_descriptions.gpt_description_generator.CellGuideConfig",
            new=MockCellGuideConfig,
        ):
            cell_counts_df = snapshot.cell_counts_cube.df[:]
            tree_builder = OntologyTreeBuilder(cell_counts_df)
            for obj in WHICH_OBJECTS_DO_NOT_EXIST:
                tree_builder.all_cell_type_ids_to_labels_in_corpus[obj["id"]] = obj["name"]
            new_gpt_descriptions = generate_new_gpt_descriptions(tree_builder.all_cell_type_ids_to_labels_in_corpus)
            self.assertTrue(compare_dicts(new_gpt_descriptions, expected_outputs))

    def test__gpt_seo_description_generator(self):
        expected_outputs = {}
        for obj in WHICH_OBJECTS_DO_NOT_EXIST:
            expected_outputs[
                obj["id"]
            ] = f"MOCK OUTPUT: {GPT_CELLTYPE_DESCRIPTION_SYSTEM_ROLE} {GPT_CELLTYPE_SEO_DESCRIPTION_USER_ROLE(obj['fake_gpt_description'])}"

        with patch(
            "backend.cellguide.pipeline.gpt_descriptions.gpt_description_generator.S3Provider", new=MockS3Provider
        ), patch(
            "backend.cellguide.pipeline.gpt_descriptions.gpt_description_generator.OpenAIProvider",
            new=MockOpenAIProvider,
        ), patch(
            "backend.cellguide.pipeline.gpt_descriptions.gpt_description_generator.CellGuideConfig",
            new=MockCellGuideConfig,
        ):
            fake_new_gpt_descriptions = {i["id"]: i["fake_gpt_description"] for i in WHICH_OBJECTS_DO_NOT_EXIST}
            new_gpt_seo_descriptions = generate_new_seo_gpt_descriptions(fake_new_gpt_descriptions)
            self.assertTrue(compare_dicts(new_gpt_seo_descriptions, expected_outputs))
