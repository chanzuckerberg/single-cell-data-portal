import json
import unittest

import pkg_resources

from dcp_prototype.backend.wrangling.migrations.common.dataset_metadata import DatasetMetadata
from dcp_prototype.backend.wrangling.migrations.common.gather_dcp_one_data import \
    generate_metadata_structure_from_targz, generate_metadata_structure_from_dir


class TestMigration(unittest.TestCase):
    def test_end_to_end_migration_from_targz(self):
        """A simple end to end test case"""
        infile = pkg_resources.resource_filename(__name__, "../../fixtures/WongAdultRetina.tar.gz")
        expectedfile = pkg_resources.resource_filename(__name__, "../../fixtures/WongAdultRetina.json")

        dataset_metadata = DatasetMetadata()
        generate_metadata_structure_from_targz(infile, dataset_metadata)
        dataset_metadata.process()
        result_project = dataset_metadata.to_dict()
        with open(expectedfile) as jfile:
            expected_project = json.load(jfile)
            title = result_project.get("projects")[0].get("title")
            self.assertEqual(title, "WongAdultRetina")
            self.assertDictEqual(result_project, expected_project)

    def test_end_to_end_migration_from_local_files(self):
        project_metadata_input_directory = pkg_resources.resource_filename(__name__, "../../fixtures/WongAdultRetina/")

        dataset_metadata = DatasetMetadata()
        generate_metadata_structure_from_dir(project_metadata_input_directory, 1, dataset_metadata)
        dataset_metadata.process()
        dictionary_representation_of_project = dataset_metadata.to_dict()

        expected_artifact_file = pkg_resources.resource_filename(__name__, "../../fixtures/WongAdultRetina.json")
        with open(expected_artifact_file) as artifact:
            expected_project = json.load(artifact)
            project_title = dictionary_representation_of_project.get("projects")[0].get("title")

            self.assertEqual(project_title, "WongAdultRetina")
            self.assertDictEqual(dictionary_representation_of_project, expected_project)
