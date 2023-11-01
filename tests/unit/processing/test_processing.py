from unittest.mock import MagicMock, patch

from backend.layers.common.entities import (
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetUploadStatus,
    DatasetValidationStatus,
)
from backend.layers.processing.process import ProcessMain
from backend.layers.processing.process_cxg import ProcessCxg
from backend.layers.processing.process_seurat import ProcessSeurat
from tests.unit.processing.base_processing_test import BaseProcessingTest, mock_config_fn


class ProcessingTest(BaseProcessingTest):
    def test_process_seurat_success(self):
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, "nothing", None, None
        )

        with patch("backend.layers.processing.process_seurat.ProcessSeurat.make_seurat") as mock:
            mock.return_value = "local.rds"
            ps = ProcessSeurat(self.business_logic, self.uri_provider, self.s3_provider)
            ps.process(dataset_version_id, "fake_bucket_name", "fake_datasets_bucket")

            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.rds_status, DatasetConversionStatus.UPLOADED)

            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.rds"))
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}.rds"))

            artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
            self.assertEqual(1, len(artifacts))
            artifact = artifacts[0]
            artifact.type = "RDS"
            artifact.uri = f"s3://fake_bucket_name/{dataset_version_id.id}/local.rds"

    def test_process_cxg_success(self):
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, "nothing", None, None
        )

        with patch("backend.layers.processing.process_cxg.ProcessCxg.make_cxg") as mock:
            mock.return_value = "local.cxg"
            ps = ProcessCxg(self.business_logic, self.uri_provider, self.s3_provider)
            ps.process(dataset_version_id, "fake_bucket_name", "fake_cxg_bucket")

            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.cxg_status, DatasetConversionStatus.UPLOADED)

            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/"))

            artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
            self.assertEqual(1, len(artifacts))
            artifact = artifacts[0]
            artifact.type = "CXG"
            artifact.uri = f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/"

    def test_reprocess_cxg_success(self):
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, "nothing", None, None
        )

        with patch("backend.layers.processing.process_cxg.ProcessCxg.make_cxg") as mock:
            mock.return_value = "local.cxg"
            ps = ProcessCxg(self.business_logic, self.uri_provider, self.s3_provider)
            self.business_logic.add_dataset_artifact(
                dataset_version_id, "h5ad", f"s3://fake_bucket_name/{dataset_id}/local.h5ad"
            )
            ps.process(dataset_version_id, "fake_bucket_name", "fake_cxg_bucket")

            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.cxg_status, DatasetConversionStatus.UPLOADED)

            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/"))

            ps.process(dataset_version_id, "fake_bucket_name", "diff_cxg_bucket", is_reprocess=True)

            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.cxg_status, DatasetConversionStatus.UPLOADED)

            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/"))

            artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
            cxg_artifact = [artifact for artifact in artifacts if artifact.type == "cxg"][0]
            self.assertTrue(cxg_artifact, f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/")

    @patch("scanpy.read_h5ad")
    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    @patch("backend.layers.processing.process_validate.ProcessValidate.extract_metadata")
    @patch("backend.layers.processing.process_seurat.ProcessSeurat.make_seurat")
    @patch("backend.layers.processing.process_cxg.ProcessCxg.make_cxg")
    def test_process_all(self, mock_cxg, mock_seurat, mock_h5ad, mock_config, mock_read_h5ad):
        mock_seurat.return_value = "local.rds"
        mock_cxg.return_value = "local.cxg"

        # Mock anndata object
        mock_anndata = MagicMock(uns=dict(), n_obs=1000, n_vars=1000)
        mock_read_h5ad.return_value = mock_anndata

        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )

        pm = ProcessMain(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        for step_name in ["download", "validate", "cxg", "seurat"]:
            assert pm.process(
                collection.version_id,
                dataset_version_id,
                step_name,
                dropbox_uri,
                "fake_bucket_name",
                "fake_datasets_bucket",
                "fake_cxg_bucket",
            )

        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/raw.h5ad"))
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"))
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}.h5ad"))
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.rds"))
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}.rds"))
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/"))

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.cxg_status, DatasetConversionStatus.UPLOADED)
        self.assertEqual(status.rds_status, DatasetConversionStatus.UPLOADED)
        self.assertEqual(status.h5ad_status, DatasetConversionStatus.UPLOADED)
        self.assertEqual(status.validation_status, DatasetValidationStatus.VALID)
        self.assertEqual(status.upload_status, DatasetUploadStatus.UPLOADED)
        self.assertEqual(status.processing_status, DatasetProcessingStatus.PENDING)
        # TODO: DatasetProcessingStatus.SUCCESS is set by a lambda that also needs to be modified. It should belong here

        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
        self.assertEqual(4, len(artifacts))
