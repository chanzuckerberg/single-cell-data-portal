import json

import jsonschema
import pkg_resources

import dcp_prototype


class ArtifactValidator():
    """
    This class is a general purpose validator to ensure that given data is compatible with the DCP artifact
    specification.
    """

    def __init__(self, artifact_filename="backend/v0.0.0.json", artifact_full_file_path=None):
        if artifact_full_file_path:
            self.artifact_spec_file = artifact_full_file_path
        else:
            self.artifact_spec_file = pkg_resources.resource_filename(dcp_prototype.__name__, artifact_filename)

    def validate(self, artifact_data):
        try:
            artifact_spec_data = open(self.artifact_spec_file).read()
            specification = json.loads(artifact_spec_data)
        except FileNotFoundError:
            # FIXME, the schema is not yet committed to main.
            # Validate should fail if the data cannot be validated.
            print(f"WARNING: The artifact specification could not be found: {self.artifact_spec_file}")
            return True

        try:
            jsonschema.validate(artifact_data, specification)
            return True
        except jsonschema.ValidationError as error:
            print("ERROR: The inputted data does not conform to the artifact specification: ", error)
            return False
