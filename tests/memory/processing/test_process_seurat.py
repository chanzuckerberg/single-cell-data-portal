"""
This script is used to test the ProcessCxg class.
"""

import shutil
import tempfile

from backend.layers.common.entities import DatasetVersionId
from backend.layers.processing.process_seurat import ProcessSeurat
from tests.unit.backend.fixtures.environment_setup import fixture_file_path

if __name__ == "__main__":
    file_name = "labeled_visium.h5ad"
    dataset_version_id = DatasetVersionId("test_dataset_id")
    with tempfile.TemporaryDirectory() as tmpdirname:
        temp_file = "/".join([tmpdirname, file_name])
        shutil.copy(fixture_file_path(file_name), temp_file)
        process = ProcessSeurat(None, None, None)
        process.make_seurat(temp_file, dataset_version_id)
