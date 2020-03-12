import unittest

import pkg_resources

from dcp_prototype.backend.wrangling.migrations.common.artifact_validation import ArtifactValidator


class TestArtifactValidator(unittest.TestCase):
    def test_validate_against_schema_successful(self):
        schema_file = pkg_resources.resource_filename(__name__, "../../fixtures/bogo_artifact_schema.json")
        sample_conformed_data = {
            "name": "John Smith",
            "nickname": "Cell Type Discoverer",
            "cell_types_discovered": 285

        }

        artifact_validator = ArtifactValidator(artifact_full_file_path=schema_file)
        valid = artifact_validator.validate(sample_conformed_data)

        self.assertTrue(valid)

    def test_validate_against_schema_fails(self):
        schema_file = pkg_resources.resource_filename(__name__, "../../fixtures/bogo_artifact_schema.json")
        sample_conformed_data = {
            "name": "John Smith",
            "nickname": "Cell Type Discoverer",
            "cell_types_discovered": "285"

        }

        artifact_validator = ArtifactValidator(artifact_full_file_path=schema_file)
        valid = artifact_validator.validate(sample_conformed_data)

        self.assertFalse(valid)
