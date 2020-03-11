import unittest
import pkg_resources
import json


from dcp_prototype.backend.wrangling.migrations.utils.migration_utils import DatasetMetadata
from dcp_prototype.backend.wrangling.migrations.utils.gather_dcp1_data import generate_metadata_structure_from_targz


class TestMigration(unittest.TestCase):
    def test_end2end(self):
        """A simple end to end test case"""
        infile = pkg_resources.resource_filename(__name__, "WongAdultRetina.tar.gz")
        expectedfile = pkg_resources.resource_filename(__name__, "WongAdultRetina.json")

        dataset_metadata = DatasetMetadata()
        generate_metadata_structure_from_targz(infile, dataset_metadata)
        dataset_metadata.process()
        result_project = dataset_metadata.to_dict()
        with open(expectedfile) as jfile:
            expected_project = json.load(jfile)
            title = result_project.get("projects")[0].get("title")
            self.assertEqual(title, "WongAdultRetina")
            self.assertDictEqual(result_project, expected_project)
