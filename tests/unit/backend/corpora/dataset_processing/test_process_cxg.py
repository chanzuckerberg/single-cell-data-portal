import pathlib
import tempfile
from unittest.mock import patch

from backend.corpora.common.corpora_orm import (
    DatasetArtifactFileType,
)

from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.dataset_processing.exceptions import ConversionFailed
from backend.corpora.dataset_processing.process_cxg import process_cxg
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestDatasetProcessing(DataPortalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tmp_dir = tempfile.mkdtemp()
        cls.h5ad_filename = pathlib.Path(cls.tmp_dir, "local.h5ad")
        cls.cxg_filename = pathlib.Path(cls.tmp_dir, "local.cxg")

        cls.h5ad_filename.touch()
        cls.cxg_filename.touch()

    @patch("backend.corpora.dataset_processing.process_cxg.make_cxg")
    @patch("backend.corpora.dataset_processing.process_cxg.subprocess.run")
    def test_create_explorer_cxg(self, mock_subprocess, mock_cxg):
        mock_cxg.return_value = str(self.cxg_filename)
        dataset = self.generate_dataset(self.session)
        dataset_id = dataset.id

        explorer_bucket = "CELLXGENE-HOSTED-TEST"

        process_cxg(str(self.h5ad_filename), dataset_id, explorer_bucket)

        dataset = Dataset.get(self.session, dataset_id)
        artifacts = dataset.artifacts

        self.assertEqual(len(artifacts), 1)
        self.assertEqual(artifacts[0].dataset_id, dataset_id)
        self.assertEqual(artifacts[0].s3_uri, f"s3://{explorer_bucket}/{artifacts[0].id}.cxg/")
        self.assertEqual(artifacts[0].filetype, DatasetArtifactFileType.CXG)

    @patch("backend.corpora.dataset_processing.process_cxg.make_cxg")
    def test_process_with_cxg_conversion_failures(self, mock_cxg):
        mock_cxg.side_effect = RuntimeError("cxg conversion failed")
        test_dataset_id = self.generate_dataset(
            self.session,
        ).id
        artifact_bucket = "test-artifact-bucket"
        with self.assertRaises(ConversionFailed):
            process_cxg(str(self.h5ad_filename), test_dataset_id, artifact_bucket)
