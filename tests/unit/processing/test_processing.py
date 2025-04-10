from unittest.mock import Mock, patch

from backend.layers.common.entities import (
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetUploadStatus,
    DatasetValidationStatus,
)
from backend.layers.common.ingestion_manifest import IngestionManifest
from backend.layers.processing.exceptions import ValidationAtacFailed
from backend.layers.processing.process import ProcessMain
from backend.layers.processing.process_cxg import ProcessCxg
from tests.unit.processing.base_processing_test import BaseProcessingTest


class ProcessingTest(BaseProcessingTest):
    def test_process_cxg_success(self):
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, "http://fake.url", None, None
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
            collection.version_id, "http://fake.url", None, None
        )

        with patch("backend.layers.processing.process_cxg.ProcessCxg.make_cxg") as mock:
            mock.return_value = "local.cxg"
            ps = ProcessCxg(self.business_logic, self.uri_provider, self.s3_provider)
            self.business_logic.add_dataset_artifact(
                dataset_version_id, DatasetArtifactType.H5AD, f"s3://fake_bucket_name/{dataset_id}/local.h5ad"
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

    @patch("backend.layers.processing.process_add_labels.ProcessAddLabels.populate_dataset_citation")
    @patch("backend.layers.processing.process_add_labels.ProcessAddLabels.extract_metadata")
    @patch("backend.layers.processing.process_cxg.ProcessCxg.make_cxg")
    def test_process_anndata(self, mock_cxg, mock_extract_h5ad, mock_dataset_citation):
        mock_cxg.return_value = "local.cxg"

        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        manifest = IngestionManifest(anndata=dropbox_uri)
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )
        self.schema_validator.check_anndata_requires_fragment.return_value = False
        pm = ProcessMain(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pm.download_from_source_uri = lambda x, y: y

        for step_name in ["validate_anndata", "validate_atac", "add_labels", "cxg"]:
            assert pm.process(
                collection.version_id,
                dataset_version_id,
                step_name,
                manifest,
                "fake_bucket_name",
                "fake_datasets_bucket",
                "fake_cxg_bucket",
            )

        dataset = self.business_logic.get_dataset_version(dataset_version_id)

        def assert_s3_matches_db(uri, artifact_type):
            self.assertTrue(self.s3_provider.uri_exists(uri))
            artifact = [artifact for artifact in dataset.artifacts if artifact.type == artifact_type]
            self.assertTrue(len(artifact) == 1)
            self.assertEqual(artifact[0].uri, uri)

        assert_s3_matches_db(f"s3://fake_bucket_name/{dataset_version_id.id}/raw.h5ad", DatasetArtifactType.RAW_H5AD)
        assert_s3_matches_db(f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad", DatasetArtifactType.H5AD)
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}.h5ad"))
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/"))
        self.assertFalse(
            self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}-fragment.tsv.bgz")
        )
        self.assertFalse(
            self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}-fragment.tsv.bgz.tbi")
        )
        self.assertFalse(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.rds"))
        self.assertFalse(self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}.rds"))

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.cxg_status, DatasetConversionStatus.UPLOADED)
        self.assertEqual(status.rds_status, DatasetConversionStatus.SKIPPED)
        self.assertEqual(status.h5ad_status, DatasetConversionStatus.UPLOADED)
        self.assertEqual(status.atac_status, DatasetConversionStatus.SKIPPED)
        self.assertEqual(status.validation_status, DatasetValidationStatus.VALID)
        self.assertEqual(status.upload_status, DatasetUploadStatus.UPLOADED)
        self.assertEqual(status.processing_status, DatasetProcessingStatus.PENDING)
        # TODO: DatasetProcessingStatus.SUCCESS is set by a lambda that also needs to be modified. It should belong here

        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
        self.assertEqual(3, len(artifacts))

    def test_process_atac_ValidationAtacFailed(self):
        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        manifest = IngestionManifest(anndata=dropbox_uri, atac_fragment=dropbox_uri)
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )
        pm = ProcessMain(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pm.process_validate_atac_seq.process = Mock(side_effect=ValidationAtacFailed(errors=["failure 1", "failure 2"]))

        assert not pm.process(
            collection.version_id,
            dataset_version_id,
            "validate_atac",
            manifest,
            "fake_bucket_name",
            "fake_datasets_bucket",
            "fake_cxg_bucket",
        )

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.validation_status, DatasetValidationStatus.INVALID)
        self.assertEqual(status.validation_message, "\n".join(["failure 1", "failure 2"]))
